% C02_batch_convert_to_bw.m
% Batch converts image sequence to black and white
% Purpose: Process all images in the dataset and save as binary images
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
output_dir = p.output_dir_02;
image_pattern = p.image_pattern;
gaussian_sigma = p.gaussian_sigma;
threshold_percentile = p.threshold_percentile;
min_dot_area = p.min_dot_area;
max_dot_area = p.max_dot_area;
save_format = p.save_format;
preview_mode = p.preview_mode;
show_progress_every = p.show_progress_every;

% ========== SETUP OUTPUT DIRECTORY ==========
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% ========== GET LIST OF IMAGES ==========
fprintf('Scanning for images...\n');
image_files = dir(fullfile(data_dir, image_pattern));
num_images = length(image_files);

if num_images == 0
    error('No images found in %s matching pattern %s', data_dir, image_pattern);
end

fprintf('Found %d images\n', num_images);

if preview_mode
    num_images = min(num_images, 10);
    fprintf('PREVIEW MODE: Processing only first %d images\n', num_images);
end

% ========== BATCH PROCESSING ==========
fprintf('\nStarting batch processing...\n');
fprintf('Parameters:\n');
fprintf('  Gaussian sigma: %d\n', gaussian_sigma);
fprintf('  Threshold percentile: %d%%\n', threshold_percentile);
fprintf('  Min dot area: %d px\n', min_dot_area);
fprintf('  Max dot area: %d px\n', max_dot_area);
fprintf('  Output format: %s\n', save_format);
fprintf('------------------------------------\n');

tic; % Start timer

for i = 1:num_images
    % Load image
    img_name = image_files(i).name;
    img_path = fullfile(data_dir, img_name);
    img = imread(img_path);

    % Convert to grayscale if needed
    if size(img, 3) == 3
        img_gray = rgb2gray(img);
    else
        img_gray = img;
    end

    % Background subtraction
    img_bg = imgaussfilt(double(img_gray), gaussian_sigma);
    img_corrected = double(img_gray) - img_bg;
    img_corrected = img_corrected - min(img_corrected(:));
    img_corrected = 255 * img_corrected / max(img_corrected(:));
    img_corrected = uint8(img_corrected);

    % Calculate threshold based on percentile (on corrected image)
    threshold_value = prctile(double(img_corrected(:)), threshold_percentile);

    % Create binary mask (1 = black dots, 0 = background)
    bw = img_corrected < threshold_value;

    % Clean up noise and small regions
    bw_cleaned = bwareaopen(bw, min_dot_area);

    % Optional: Filter by maximum area
    if max_dot_area < inf
        labeled = bwlabel(bw_cleaned);
        stats = regionprops(labeled, 'Area');
        valid_labels = find([stats.Area] <= max_dot_area);
        bw_final = ismember(labeled, valid_labels);
    else
        bw_final = bw_cleaned;
    end

    % Save processed image
    % Format: prefix_bw_000001.png (index at end for video generation)
    output_name = sprintf('%s_bw_%06d.%s', p.dataset_name, i, save_format);
    output_path = fullfile(output_dir, output_name);
    imwrite(bw_final, output_path);

    % Show progress
    if mod(i, show_progress_every) == 0 || i == 1 || i == num_images
        elapsed = toc;
        avg_time = elapsed / i;
        remaining = avg_time * (num_images - i);
        fprintf('[%4d/%4d] Processed: %s (%.2f s/frame, ETA: %.1f s)\n', ...
            i, num_images, img_name, avg_time, remaining);
    end
end

total_time = toc;

% ========== SUMMARY ==========
fprintf('\n========== PROCESSING COMPLETE ==========\n');
fprintf('Total images processed: %d\n', num_images);
fprintf('Total time: %.2f seconds\n', total_time);
fprintf('Average time per image: %.3f seconds\n', total_time / num_images);
fprintf('Output directory: %s\n', output_dir);
fprintf('=========================================\n');

% ========== SHOW SAMPLE RESULTS ==========
fprintf('\nGenerating sample comparison...\n');

% Show first, middle, and last frames
sample_indices = [1, floor(num_images/2), num_images];
figure('Name', 'C02 - Sample Results', 'Position', [100, 100, 1400, 600]);

for idx = 1:3
    i = sample_indices(idx);

    % Load original
    img_name = image_files(i).name;
    img_path = fullfile(data_dir, img_name);
    img = imread(img_path);

    % Load processed
    output_name = sprintf('%s_bw_%06d.%s', p.dataset_name, i, save_format);
    output_path = fullfile(output_dir, output_name);
    bw = imread(output_path);

    % Show original
    subplot(2, 3, idx);
    imshow(img);
    title(sprintf('Original - Frame %d', i));

    % Show processed
    subplot(2, 3, idx + 3);
    imshow(bw);
    title(sprintf('Binary - Frame %d', i));
end

% ========== GENERATE VIDEO ==========
fprintf('\nGenerating video from image sequence...\n');

video_filename = fullfile(output_dir, sprintf('%s_bw_video.mp4', p.dataset_name));
video_writer = VideoWriter(video_filename, 'MPEG-4');
video_writer.FrameRate = 30;
open(video_writer);

for i = 1:num_images
    output_name = sprintf('%s_bw_%06d.%s', p.dataset_name, i, save_format);
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

    if mod(i, 50) == 0 || i == num_images
        fprintf('  Encoded %d/%d frames\n', i, num_images);
    end
end

close(video_writer);
fprintf('Video saved: %s\n', video_filename);

fprintf('\nC02 complete! Check the output directory for results.\n');
fprintf('Review the sample comparison figure to verify quality.\n');
