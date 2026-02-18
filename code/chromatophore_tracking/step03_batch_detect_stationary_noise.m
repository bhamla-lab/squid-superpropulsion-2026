% C03_detect_stationary_noise.m
% Detects stationary light noise and subtracts it from processed images
% Purpose: Remove stationary artifacts/noise that appear in most frames
% Author: Claude
% Date: 2025-10-11

clear; close all; clc;

% ========== SELECT DATASET ==========
% Change this number to select different dataset:
%   1 = first dataset, 2 = second dataset, etc.
DATASET_INDEX = 9;

% ========== LOAD SHARED PARAMETERS ==========
p = C00_parameters(DATASET_INDEX);

% Extract parameters
input_dir = p.output_dir_02;  % Use processed images from C02
image_pattern = sprintf('*.%s', p.save_format);

% Output directory for cleaned images
output_dir = p.output_dir_03;

% ========== STATIONARY DETECTION PARAMETERS ==========
num_frames_to_check = p.num_frames_to_check;
stationary_threshold = p.stationary_threshold;
flicker_threshold = p.flicker_threshold;
flicker_max_size = p.flicker_max_size;

% ========== SETUP OUTPUT DIRECTORY ==========
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% ========== LOAD PROCESSED IMAGES ==========
fprintf('Loading processed images...\n');
image_files = dir(fullfile(input_dir, image_pattern));
total_images = length(image_files);

if total_images == 0
    error('No processed images found in %s\n', input_dir);
end

fprintf('Found %d processed images\n', total_images);

% Determine number of frames for noise detection
if num_frames_to_check > 0
    num_images = min(total_images, num_frames_to_check);
else
    num_images = total_images;
end

fprintf('Analyzing first %d frames for stationary noise...\n', num_images);

% ========== ACCUMULATE DOTS ACROSS FRAMES ==========
fprintf('Accumulating dots across frames...\n');

% Load first image to get dimensions
first_img = imread(fullfile(input_dir, image_files(1).name));
[height, width] = size(first_img);

% Initialize accumulator (counts how many times each pixel is a dot)
accumulator = zeros(height, width, 'uint16');

% Accumulate all frames
tic;
for i = 1:num_images
    img_path = fullfile(input_dir, image_files(i).name);
    img = imread(img_path);

    % Add to accumulator (binary image: 0 or 1)
    accumulator = accumulator + uint16(img > 0);

    if mod(i, 50) == 0 || i == num_images
        fprintf('  Processed %d/%d frames\n', i, num_images);
    end
end
elapsed = toc;

fprintf('Accumulation complete in %.2f seconds\n', elapsed);

% ========== DETECT STATIONARY DOTS ==========
fprintf('\nDetecting stationary dots...\n');

% Calculate appearance frequency for each pixel
appearance_frequency = double(accumulator) / num_images;

% Identify stationary dots (appear in >= threshold fraction of frames)
stationary_mask_initial = appearance_frequency >= stationary_threshold;

% Count initial stationary pixels
num_stationary_pixels_initial = sum(stationary_mask_initial(:));

fprintf('Stationary pixels (>= %.0f%% frames): %d\n', ...
    stationary_threshold * 100, num_stationary_pixels_initial);

% ========== EXPAND TO FULL CONNECTED REGIONS ==========
fprintf('Expanding to remove entire connected components...\n');

% For each frame, find connected components that overlap with stationary pixels
% Then mark ALL pixels in those components as stationary
stationary_mask_expanded = false(size(stationary_mask_initial));

for i = 1:num_images
    % Load frame
    img_path = fullfile(input_dir, image_files(i).name);
    img = imread(img_path);

    % Find connected components in this frame
    cc = bwconncomp(img, 4);

    % Check each component
    for j = 1:cc.NumObjects
        component_pixels = cc.PixelIdxList{j};

        % If ANY pixel in this component is stationary, mark ALL pixels as stationary
        if any(stationary_mask_initial(component_pixels))
            stationary_mask_expanded(component_pixels) = true;
        end
    end
end

% Use expanded mask
stationary_mask = stationary_mask_expanded;
num_stationary_pixels = sum(stationary_mask(:));

fprintf('Expanded stationary pixels: %d (removed %d additional partial regions)\n', ...
    num_stationary_pixels, num_stationary_pixels - num_stationary_pixels_initial);

% ========== FILTER FLICKERING SMALL DOTS ==========
fprintf('Detecting flickering small dots...\n');

% Detect pixels with frequency > flicker_threshold and < stationary_threshold
flicker_mask = (appearance_frequency > flicker_threshold) & (appearance_frequency < stationary_threshold);

% Find connected components of flickering pixels
cc_flicker = bwconncomp(flicker_mask, 4);
flicker_small_mask = false(size(flicker_mask));

% Keep only components with size <= flicker_max_size
num_flicker_removed = 0;
for j = 1:cc_flicker.NumObjects
    component_pixels = cc_flicker.PixelIdxList{j};
    component_size = length(component_pixels);

    if component_size <= flicker_max_size
        flicker_small_mask(component_pixels) = true;
        num_flicker_removed = num_flicker_removed + component_size;
    end
end

% Add flickering small dots to stationary mask
stationary_mask = stationary_mask | flicker_small_mask;
num_stationary_pixels = sum(stationary_mask(:));

fprintf('Flickering small dots removed: %d pixels\n', num_flicker_removed);

% Label connected components in final mask
[labeled, num_stationary_regions] = bwlabel(stationary_mask);
stats = regionprops(labeled, 'Area', 'Centroid', 'BoundingBox');

fprintf('Stationary regions: %d\n', num_stationary_regions);

% ========== SUBTRACT STATIONARY NOISE FROM ALL FRAMES ==========
fprintf('\nSubtracting stationary noise from all frames...\n');

tic;
for i = 1:total_images
    % Load original processed image
    img_path = fullfile(input_dir, image_files(i).name);
    img = imread(img_path);

    % Subtract expanded stationary mask
    img_cleaned = img & ~stationary_mask;

    % Save cleaned image
    % Format: prefix_cleaned_000001.png (index at end for video generation)
    output_name = sprintf('%s_cleaned_%06d.%s', p.dataset_name, i, p.save_format);
    output_path = fullfile(output_dir, output_name);
    imwrite(img_cleaned, output_path);

    % Show progress
    if mod(i, 50) == 0 || i == 1 || i == total_images
        fprintf('  Cleaned %d/%d frames\n', i, total_images);
    end
end
elapsed_clean = toc;

fprintf('Noise subtraction complete in %.2f seconds\n', elapsed_clean);

% ========== VISUALIZE RESULTS ==========
fprintf('\nGenerating visualization...\n');

% Load cleaned first frame for comparison
first_cleaned_path = fullfile(output_dir, sprintf('%s_cleaned_%06d.%s', ...
    p.dataset_name, 1, p.save_format));
first_cleaned = imread(first_cleaned_path);

figure('Name', 'C03 - Stationary Noise Subtraction', 'Position', [100, 100, 1600, 900]);

% Show first frame (before cleaning)
subplot(2, 4, 1);
imshow(first_img);
title('1. Before Cleaning');

% Show stationary mask
subplot(2, 4, 2);
imshow(stationary_mask);
title(sprintf('2. Stationary Mask (N=%d)', num_stationary_pixels));

% Show first frame (after cleaning)
subplot(2, 4, 3);
imshow(first_cleaned);
title('3. After Cleaning');

% Show what was removed
subplot(2, 4, 4);
removed = first_img & stationary_mask;
imshow(removed);
title('4. Removed Noise');

% Show accumulator (heat map)
subplot(2, 4, 5);
imagesc(accumulator);
colorbar;
colormap('hot');
title(sprintf('Accumulator (Max=%d)', max(accumulator(:))));
axis image;

% Show appearance frequency
subplot(2, 4, 6);
imagesc(appearance_frequency);
colorbar;
colormap('hot');
caxis([0 1]);
title('Appearance Frequency');
axis image;

% Show overlay on first frame
subplot(2, 4, 7);
imshow(first_img);
hold on;
[y, x] = find(stationary_mask);
plot(x, y, 'r.', 'MarkerSize', 2);
title('Stationary Dots Overlay');

% Show stationary regions with labels
subplot(2, 4, 8);
imshow(first_img);
hold on;
for i = 1:min(length(stats), 20)  % Limit to 20 labels for clarity
    bbox = stats(i).BoundingBox;
    rectangle('Position', bbox, 'EdgeColor', 'r', 'LineWidth', 1.5);
    centroid = stats(i).Centroid;
    plot(centroid(1), centroid(2), 'g+', 'MarkerSize', 8, 'LineWidth', 2);
end
title(sprintf('Stationary Regions (N=%d)', num_stationary_regions));

% ========== PRINT STATISTICS ==========
fprintf('\n========== SUMMARY ==========\n');
fprintf('Total frames processed: %d\n', total_images);
fprintf('Frames analyzed for noise: %d\n', num_images);
fprintf('Stationary threshold: %.0f%%\n', stationary_threshold * 100);
fprintf('Stationary pixels removed: %d\n', num_stationary_pixels);
fprintf('Stationary regions: %d\n', num_stationary_regions);
fprintf('Output directory: %s\n', output_dir);

if num_stationary_regions > 0
    areas = [stats.Area];
    fprintf('\nStationary region statistics:\n');
    fprintf('  Min area: %d px\n', min(areas));
    fprintf('  Max area: %d px\n', max(areas));
    fprintf('  Mean area: %.1f px\n', mean(areas));
    fprintf('  Median area: %.1f px\n', median(areas));

    fprintf('\nTop 10 largest stationary regions:\n');
    [sorted_areas, idx] = sort(areas, 'descend');
    for i = 1:min(10, length(sorted_areas))
        region_id = idx(i);
        fprintf('  Region %d: %d px at (%.1f, %.1f)\n', ...
            region_id, sorted_areas(i), ...
            stats(region_id).Centroid(1), stats(region_id).Centroid(2));
    end
end
fprintf('=============================\n');

% ========== GENERATE VIDEO ==========
fprintf('\nGenerating video from cleaned image sequence...\n');

video_filename = fullfile(output_dir, sprintf('%s_cleaned_video.mp4', p.dataset_name));
video_writer = VideoWriter(video_filename, 'MPEG-4');
video_writer.FrameRate = 30;
open(video_writer);

for i = 1:total_images
    output_name = sprintf('%s_cleaned_%06d.%s', p.dataset_name, i, p.save_format);
    output_path = fullfile(output_dir, output_name);
    frame = imread(output_path);

    % Convert logical/binary to uint8
    if islogical(frame)
        frame = uint8(frame) * 255;
    end

    % Convert grayscale to RGB if needed
    if size(frame, 3) == 1
        frame = repmat(frame, [1, 1, 3]);
    end

    writeVideo(video_writer, frame);

    if mod(i, 50) == 0 || i == total_images
        fprintf('  Encoded %d/%d frames\n', i, total_images);
    end
end

close(video_writer);
fprintf('Video saved: %s\n', video_filename);

fprintf('\nC03 complete!\n');
fprintf('- Stationary noise detected and removed from all %d frames\n', total_images);
fprintf('- Cleaned images saved to: %s\n', output_dir);
fprintf('- Review the figure to verify noise removal\n');
