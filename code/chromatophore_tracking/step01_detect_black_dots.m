% C01_detect_black_dots.m
% Detects black dots in the first image of the sequence
% Purpose: Test detection before batch processing
% Author: Claude
% Date: 2025-10-11

clear; close all; clc;

% ========== SELECT DATASET ==========
% Change this number to select different dataset:
%   1 = first dataset, 2 = second dataset, etc.
DATASET_INDEX = 9;

% ========== LOAD SHARED PARAMETERS ==========
p = C00_parameters(DATASET_INDEX);

% Extract parameters for easier access
data_dir = p.data_dir;
first_image_name = p.first_image_name;
gaussian_sigma = p.gaussian_sigma;
threshold_percentile = p.threshold_percentile;
min_dot_area = p.min_dot_area;
max_dot_area = p.max_dot_area;

% ========== LOAD FIRST IMAGE ==========
fprintf('Loading first image...\n');
img_path = fullfile(data_dir, first_image_name);

if ~exist(img_path, 'file')
    error('Image file not found: %s', img_path);
end

img = imread(img_path);
fprintf('Image loaded: %dx%d\n', size(img, 1), size(img, 2));

% Convert to grayscale if needed
if size(img, 3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end

% ========== BACKGROUND SUBTRACTION ==========
fprintf('Subtracting background...\n');

% Create Gaussian-blurred background
img_bg = imgaussfilt(double(img_gray), gaussian_sigma);

% Subtract background
img_corrected = double(img_gray) - img_bg;

% Normalize to 0-255 range
img_corrected = img_corrected - min(img_corrected(:));
img_corrected = 255 * img_corrected / max(img_corrected(:));
img_corrected = uint8(img_corrected);

fprintf('Background subtraction complete (sigma = %d)\n', gaussian_sigma);

% ========== DETECT BLACK DOTS ==========
fprintf('Detecting black dots...\n');

% Calculate threshold based on percentile (on corrected image)
threshold_value = prctile(double(img_corrected(:)), threshold_percentile);
fprintf('Threshold value (%.1f percentile): %.1f\n', threshold_percentile, threshold_value);

% Create binary mask (1 = black dots, 0 = background)
bw = img_corrected < threshold_value;

% Clean up noise and small regions
bw_cleaned = bwareaopen(bw, min_dot_area);

% Label connected components
[labeled, num_dots] = bwlabel(bw_cleaned);

% Filter by area (both min and max)
stats = regionprops(labeled, 'Area', 'Centroid', 'BoundingBox');
valid_dots = [stats.Area] >= min_dot_area & [stats.Area] <= max_dot_area;
num_valid_dots = sum(valid_dots);

% Create final binary image with only valid dots
valid_labels = find(valid_dots);
bw_final = ismember(labeled, valid_labels);

fprintf('Total regions after min filter: %d\n', num_dots);
fprintf('Valid dots (area %d-%d px): %d\n', min_dot_area, max_dot_area, num_valid_dots);

% ========== VISUALIZE RESULTS ==========
figure('Name', 'C01 - Black Dot Detection', 'Position', [100, 100, 1400, 900]);

% Original image
subplot(3, 3, 1);
imshow(img);
title('1. Original Image');

% Grayscale
subplot(3, 3, 2);
imshow(img_gray);
title('2. Grayscale');

% Background (blurred)
subplot(3, 3, 3);
imshow(uint8(img_bg));
title(sprintf('3. Background (σ=%d)', gaussian_sigma));

% Background-corrected image
subplot(3, 3, 4);
imshow(img_corrected);
title('4. Background Subtracted');

% Histogram comparison
subplot(3, 3, 5);
histogram(double(img_gray(:)), 100, 'FaceAlpha', 0.5, 'DisplayName', 'Original');
hold on;
histogram(double(img_corrected(:)), 100, 'FaceAlpha', 0.5, 'DisplayName', 'Corrected');
xlabel('Intensity');
ylabel('Frequency');
title('5. Histogram Comparison');
legend('Location', 'best');
grid on;

% Histogram with threshold
subplot(3, 3, 6);
histogram(double(img_corrected(:)), 100);
hold on;
xline(threshold_value, 'r', 'LineWidth', 2);
xlabel('Intensity');
ylabel('Frequency');
title(sprintf('6. Threshold = %.1f (%d%%)', threshold_value, threshold_percentile));
legend('Corrected Image', 'Threshold');
grid on;

% Binary mask (before cleanup)
subplot(3, 3, 7);
imshow(bw);
title('7. Binary (Before Cleanup)');

% Binary mask (after min area filter)
subplot(3, 3, 8);
imshow(bw_cleaned);
title(sprintf('8. Min Area Filter (N=%d)', num_dots));

% Final result (min + max area filter)
subplot(3, 3, 9);
imshow(bw_final);
title(sprintf('9. Final Result (N=%d)', num_valid_dots));

% Print summary statistics
fprintf('\n========== SUMMARY ==========\n');
if num_valid_dots > 0
    areas = [stats(valid_dots).Area];
    fprintf('Dot area statistics:\n');
    fprintf('  Min area: %.1f px\n', min(areas));
    fprintf('  Max area: %.1f px\n', max(areas));
    fprintf('  Mean area: %.1f px\n', mean(areas));
    fprintf('  Median area: %.1f px\n', median(areas));
end
fprintf('=============================\n');

fprintf('\nC01 complete! Review the figure to verify detection.\n');
fprintf('If results look good, proceed to C02 for batch processing.\n');
