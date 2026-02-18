% C04_track_dots.m
% Tracks each dot across the image sequence
% Purpose: Link dots between frames to create trajectories
% Author: Claude
% Date: 2025-10-11

clear; close all; clc;

% ========== SELECT DATASET ==========
% Change this number to select different dataset:
%   1 = first dataset, 2 = second dataset, etc.
DATASET_INDEX = 9;

% ========== LOAD SHARED PARAMETERS ==========
p = C00_parameters(DATASET_INDEX);

% Input/output directories
input_dir = p.output_dir_03;
output_dir = p.output_dir_04;
image_pattern = sprintf('*.%s', p.save_format);

% ========== TRACKING PARAMETERS (from C00) ==========
max_displacement = p.max_displacement;
min_track_length = p.min_track_length_c04;
num_frames_to_process = p.num_frames_to_process;

% ========== SETUP OUTPUT DIRECTORY ==========
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% ========== LOAD IMAGES ==========
fprintf('Loading cleaned images...\n');
image_files = dir(fullfile(input_dir, image_pattern));
total_images = length(image_files);

if total_images == 0
    error('No cleaned images found in %s\n', input_dir);
end

if num_frames_to_process > 0
    total_images = min(total_images, num_frames_to_process);
end

fprintf('Found %d cleaned images to track\n', total_images);

% ========== DETECT DOTS IN ALL FRAMES ==========
fprintf('\nDetecting dots in all frames...\n');

% Initialize storage for all dots
all_dots = cell(total_images, 1);

tic;
for i = 1:total_images
    % Load image
    img_path = fullfile(input_dir, image_files(i).name);
    img = imread(img_path);

    % Find connected components (4-connectivity: corners don't connect)
    cc = bwconncomp(img, 4);
    stats = regionprops(cc, 'Centroid', 'Area');

    % Store centroids and areas
    if ~isempty(stats)
        centroids = vertcat(stats.Centroid);
        areas = [stats.Area]';
        all_dots{i} = [centroids, areas];  % [x, y, area]
    else
        all_dots{i} = zeros(0, 3);
    end

    if mod(i, 50) == 0 || i == total_images
        fprintf('  Detected dots in frame %d/%d\n', i, total_images);
    end
end
elapsed = toc;

fprintf('Dot detection complete in %.2f seconds\n', elapsed);

% Print dot count summary
num_dots_per_frame = cellfun(@(x) size(x, 1), all_dots);
fprintf('Dots per frame: min=%d, max=%d, mean=%.1f\n', ...
    min(num_dots_per_frame), max(num_dots_per_frame), mean(num_dots_per_frame));

% ========== TRACK DOTS ACROSS FRAMES ==========
fprintf('\nTracking dots across frames...\n');

% Initialize tracking variables
tracks = {};  % Cell array of tracks, each track is [frame, x, y, area]
active_tracks = [];  % Currently active tracks [last_x, last_y, track_id, last_frame]

track_id_counter = 0;

tic;
for frame = 1:total_images
    current_dots = all_dots{frame};
    num_dots = size(current_dots, 1);

    % Get currently active tracks
    num_active = size(active_tracks, 1);

    if num_active == 0 && num_dots > 0
        % First frame or no active tracks: start new tracks for all dots
        for d = 1:num_dots
            track_id_counter = track_id_counter + 1;
            tracks{track_id_counter} = [frame, current_dots(d, :)];
            active_tracks = [active_tracks; current_dots(d, 1:2), track_id_counter, frame];
        end
    elseif num_active > 0 && num_dots > 0
        % Match current dots to active tracks using nearest neighbor
        matched = false(num_dots, 1);
        updated_tracks = [];

        for a = 1:num_active
            last_pos = active_tracks(a, 1:2);
            track_id = active_tracks(a, 3);

            % Find nearest dot
            distances = sqrt(sum((current_dots(:, 1:2) - last_pos).^2, 2));
            [min_dist, min_idx] = min(distances);

            if min_dist <= max_displacement && ~matched(min_idx)
                % Match found: update track
                tracks{track_id} = [tracks{track_id}; frame, current_dots(min_idx, :)];
                updated_tracks = [updated_tracks; current_dots(min_idx, 1:2), track_id, frame];
                matched(min_idx) = true;
            end
        end

        % Start new tracks for unmatched dots
        for d = 1:num_dots
            if ~matched(d)
                track_id_counter = track_id_counter + 1;
                tracks{track_id_counter} = [frame, current_dots(d, :)];
                updated_tracks = [updated_tracks; current_dots(d, 1:2), track_id_counter, frame];
            end
        end

        active_tracks = updated_tracks;
    end

    if mod(frame, 50) == 0 || frame == total_images
        fprintf('  Tracked frame %d/%d (active tracks: %d, total tracks: %d)\n', ...
            frame, total_images, size(active_tracks, 1), track_id_counter);
    end
end
elapsed_track = toc;

fprintf('Tracking complete in %.2f seconds\n', elapsed_track);
fprintf('Total tracks created: %d\n', track_id_counter);

% ========== FILTER SHORT TRACKS ==========
fprintf('\nFiltering short tracks (min length = %d frames)...\n', min_track_length);

track_lengths = cellfun(@(x) size(x, 1), tracks);
valid_tracks = track_lengths >= min_track_length;
tracks_filtered = tracks(valid_tracks);
num_valid_tracks = sum(valid_tracks);

fprintf('Valid tracks: %d (%.1f%%)\n', num_valid_tracks, 100 * num_valid_tracks / track_id_counter);

% ========== SAVE TRACKING RESULTS ==========
fprintf('\nSaving tracking results...\n');

% Save all tracks to MAT file
save(fullfile(output_dir, 'tracks.mat'), 'tracks_filtered', 'track_lengths', ...
    'all_dots', 'max_displacement', 'min_track_length', 'total_images');

fprintf('Saved tracks to: %s\n', fullfile(output_dir, 'tracks.mat'));

% ========== STATISTICS ==========
fprintf('\n========== SUMMARY ==========\n');
fprintf('Frames processed: %d\n', total_images);
fprintf('Max displacement: %d pixels\n', max_displacement);
fprintf('Min track length: %d frames\n', min_track_length);
fprintf('Total tracks: %d\n', track_id_counter);
fprintf('Valid tracks: %d\n', num_valid_tracks);

if num_valid_tracks > 0
    valid_lengths = track_lengths(valid_tracks);
    fprintf('\nTrack length statistics:\n');
    fprintf('  Min: %d frames\n', min(valid_lengths));
    fprintf('  Max: %d frames\n', max(valid_lengths));
    fprintf('  Mean: %.1f frames\n', mean(valid_lengths));
    fprintf('  Median: %.1f frames\n', median(valid_lengths));
end
fprintf('=============================\n');

% ========== VISUALIZE SAMPLE TRACKS ==========
fprintf('\nGenerating visualization...\n');

% Load first and last frames
first_img = imread(fullfile(input_dir, image_files(1).name));
last_img = imread(fullfile(input_dir, image_files(total_images).name));

figure('Name', 'C04 - Dot Tracking', 'Position', [100, 100, 1600, 800]);

% Show first frame with all starting positions
subplot(1, 3, 1);
imshow(first_img);
hold on;
for i = 1:num_valid_tracks
    track = tracks_filtered{i};
    if track(1, 1) == 1  % Started in first frame
        plot(track(1, 2), track(1, 3), 'go', 'MarkerSize', 6, 'LineWidth', 1.5);
    end
end
title(sprintf('First Frame (N=%d dots)', size(all_dots{1}, 1)));

% Show all tracks overlaid
subplot(1, 3, 2);
imshow(first_img);
hold on;
colors = jet(min(num_valid_tracks, 100));  % Limit colors
for i = 1:min(num_valid_tracks, 100)  % Show first 100 tracks
    track = tracks_filtered{i};
    plot(track(:, 2), track(:, 3), '-', 'Color', colors(i, :), 'LineWidth', 1);
    plot(track(1, 2), track(1, 3), 'o', 'Color', colors(i, :), 'MarkerSize', 4);
end
title(sprintf('Trajectories (showing %d tracks)', min(num_valid_tracks, 100)));

% Show last frame with all ending positions
subplot(1, 3, 3);
imshow(last_img);
hold on;
for i = 1:num_valid_tracks
    track = tracks_filtered{i};
    if track(end, 1) == total_images  % Ended in last frame
        plot(track(end, 2), track(end, 3), 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    end
end
title(sprintf('Last Frame (N=%d dots)', size(all_dots{total_images}, 1)));

fprintf('\nC04 complete!\n');
fprintf('- Tracked %d valid trajectories\n', num_valid_tracks);
fprintf('- Results saved to: %s\n', output_dir);
fprintf('- Review the figure to verify tracking quality\n');
