% C05_review_tracking.m
% Reviews tracking results with interactive polygon region selection
% Purpose: Select region and show only tracks starting in that region
% Author: Claude
% Date: 2025-10-11

clear; close all; clc;

% ========== SELECT DATASET ==========
% Change this number to select different dataset:
%   1 = first dataset, 2 = second dataset, etc.
DATASET_INDEX = 9;

% ========== LOAD SHARED PARAMETERS ==========
p = C00_parameters(DATASET_INDEX);

% ========== CONFIGURATION ==========
tracking_dir = p.output_dir_04;
tracks_file = fullfile(tracking_dir, 'tracks.mat');
polygon_file = fullfile(tracking_dir, 'roi_polygon.mat');
image_dir = p.output_dir_03;
image_pattern = sprintf('*.%s', p.save_format);

% Filtering parameters (from C00)
max_displacement_per_frame = p.max_displacement_per_frame;
min_track_length_filter = p.min_track_length_c05;
min_total_displacement = p.min_total_displacement;

% ========== LOAD TRACKING DATA ==========
fprintf('Loading tracking data...\n');

if ~exist(tracks_file, 'file')
    error('Tracking data not found: %s\nPlease run C04 first.', tracks_file);
end

load(tracks_file, 'tracks_filtered');
num_tracks = length(tracks_filtered);
fprintf('Loaded %d tracks\n', num_tracks);

% Load first frame
image_files = dir(fullfile(image_dir, image_pattern));
first_img = imread(fullfile(image_dir, image_files(1).name));

% ========== POLYGON REGION SELECTION ==========
fprintf('\n========== REGION SELECTION ==========\n');

% Check for saved polygons
polygon_dir = tracking_dir;
saved_polygons = dir(fullfile(polygon_dir, 'roi_polygon_*.mat'));

if ~isempty(saved_polygons)
    fprintf('Found %d saved polygon(s):\n', length(saved_polygons));
    for i = 1:length(saved_polygons)
        fprintf('  [%d] %s\n', i, saved_polygons(i).name);
    end
    fprintf('  [%d] Draw new polygon\n', length(saved_polygons) + 1);

    % Ask user to select
    user_choice = input(sprintf('Select polygon (1-%d): ', length(saved_polygons) + 1));

    if user_choice >= 1 && user_choice <= length(saved_polygons)
        % Load selected polygon
        selected_file = fullfile(polygon_dir, saved_polygons(user_choice).name);
        load(selected_file, 'roi_x', 'roi_y', 'polygon_name');
        if ~exist('polygon_name', 'var')
            polygon_name = saved_polygons(user_choice).name;
        end
        fprintf('Loaded polygon: %s\n', polygon_name);
    else
        % Draw new polygon
        polygon_name = input('Enter name for new polygon: ', 's');
        if isempty(polygon_name)
            polygon_name = sprintf('polygon_%s', datestr(now, 'yyyymmdd_HHMMSS'));
        end

        fprintf('Draw polygon on image. Double-click to finish.\n');
        figure('Name', sprintf('Select ROI Polygon: %s', polygon_name));
        imshow(first_img);
        title(sprintf('Draw polygon: %s (Double-click to finish)', polygon_name));
        h = drawpolygon();
        roi_x = h.Position(:, 1);
        roi_y = h.Position(:, 2);
        close(gcf);

        % Save with name
        polygon_filename = sprintf('roi_polygon_%s.mat', strrep(polygon_name, ' ', '_'));
        polygon_filepath = fullfile(polygon_dir, polygon_filename);
        save(polygon_filepath, 'roi_x', 'roi_y', 'polygon_name');
        fprintf('Saved new polygon: %s\n', polygon_filename);
    end
else
    % No saved polygons, draw new one
    fprintf('No saved polygons found.\n');
    polygon_name = input('Enter name for new polygon: ', 's');
    if isempty(polygon_name)
        polygon_name = sprintf('polygon_%s', datestr(now, 'yyyymmdd_HHMMSS'));
    end

    fprintf('Draw polygon on image. Double-click to finish.\n');
    figure('Name', sprintf('Select ROI Polygon: %s', polygon_name));
    imshow(first_img);
    title(sprintf('Draw polygon: %s (Double-click to finish)', polygon_name));
    h = drawpolygon();
    roi_x = h.Position(:, 1);
    roi_y = h.Position(:, 2);
    close(gcf);

    % Save with name
    polygon_filename = sprintf('roi_polygon_%s.mat', strrep(polygon_name, ' ', '_'));
    polygon_filepath = fullfile(polygon_dir, polygon_filename);
    save(polygon_filepath, 'roi_x', 'roi_y', 'polygon_name');
    fprintf('Saved polygon: %s\n', polygon_filename);
end

% ========== FILTER TRACKS BY POLYGON ==========
fprintf('\nFiltering tracks by polygon...\n');

tracks_in_roi = false(num_tracks, 1);

for i = 1:num_tracks
    track = tracks_filtered{i};
    start_x = track(1, 2);
    start_y = track(1, 3);

    % Check if starting position is inside polygon
    tracks_in_roi(i) = inpolygon(start_x, start_y, roi_x, roi_y);
end

tracks_roi = tracks_filtered(tracks_in_roi);
num_tracks_roi = sum(tracks_in_roi);

fprintf('Tracks in ROI: %d / %d (%.1f%%)\n', ...
    num_tracks_roi, num_tracks, 100 * num_tracks_roi / num_tracks);

% ========== FILTER BY CONTINUITY (SPATIAL JUMPS ONLY) ==========
fprintf('\nFiltering by spatial continuity (removing tracks with large jumps)...\n');

continuous_tracks_mask = true(num_tracks_roi, 1);

for i = 1:num_tracks_roi
    track = tracks_roi{i};

    % Check for spatial jumps (large displacement between consecutive positions)
    positions = track(:, 2:3);
    if size(positions, 1) > 1
        % Calculate displacement between all consecutive positions in track
        displacements = sqrt(sum(diff(positions).^2, 2));

        % If any displacement exceeds threshold, reject track
        if any(displacements > max_displacement_per_frame)
            continuous_tracks_mask(i) = false;
        end
    end
end

tracks_roi_continuous = tracks_roi(continuous_tracks_mask);
tracks_roi_discontinuous = tracks_roi(~continuous_tracks_mask);

num_continuous = sum(continuous_tracks_mask);
num_discontinuous = sum(~continuous_tracks_mask);

fprintf('Spatially continuous tracks: %d\n', num_continuous);
fprintf('Tracks with jumps: %d (rejected)\n', num_discontinuous);

% ========== FILTER BY TRACK LENGTH ==========
fprintf('\nFiltering by track length...\n');

track_lengths_continuous = cellfun(@(x) size(x, 1), tracks_roi_continuous);
long_tracks_mask = track_lengths_continuous >= min_track_length_filter;

tracks_roi_short = tracks_roi_continuous(~long_tracks_mask);
tracks_roi_long = tracks_roi_continuous(long_tracks_mask);

num_short = sum(~long_tracks_mask);
num_long = sum(long_tracks_mask);

fprintf('Short continuous tracks (< %d frames): %d (rejected, shown in gray)\n', ...
    min_track_length_filter, num_short);
fprintf('Long continuous tracks (>= %d frames): %d (kept for analysis)\n', ...
    min_track_length_filter, num_long);

% ========== FILTER STATIONARY DOTS ==========
fprintf('\nFiltering stationary dots (dots that barely move)...\n');

moving_tracks_mask = true(num_long, 1);

for i = 1:num_long
    track = tracks_roi_long{i};
    positions = track(:, 2:3);

    % Calculate total displacement from start to end
    start_pos = positions(1, :);
    end_pos = positions(end, :);
    total_displacement = sqrt(sum((end_pos - start_pos).^2));

    % Reject if total displacement is too small
    if total_displacement < min_total_displacement
        moving_tracks_mask(i) = false;
    end
end

tracks_roi_stationary = tracks_roi_long(~moving_tracks_mask);
tracks_roi_moving = tracks_roi_long(moving_tracks_mask);

num_stationary = sum(~moving_tracks_mask);
num_moving = sum(moving_tracks_mask);

fprintf('Stationary tracks (< %d pixels total displacement): %d (rejected)\n', ...
    min_total_displacement, num_stationary);
fprintf('Moving tracks (>= %d pixels total displacement): %d (KEPT)\n', ...
    min_total_displacement, num_moving);

% Update for visualization - use moving tracks as the final filtered set
tracks_roi_long = tracks_roi_moving;
num_long = num_moving;

% ========== SAVE RESULTS ==========
fprintf('\nSaving results to %s...\n', tracking_dir);

% Create output filename with polygon name
output_filename = sprintf('filtered_tracks_%s.mat', strrep(polygon_name, ' ', '_'));
output_file = fullfile(tracking_dir, output_filename);

save(output_file, 'tracks_roi_long', 'tracks_roi_short', 'roi_x', 'roi_y', ...
    'polygon_name', 'min_track_length_filter', 'num_long', 'num_short');

fprintf('Saved to: %s\n', output_file);

% ========== PRINT STATISTICS ==========
fprintf('\n========== SUMMARY ==========\n');
fprintf('Polygon: %s\n', polygon_name);
fprintf('Total tracks: %d\n', num_tracks);
fprintf('Tracks in ROI: %d (%.1f%%)\n', num_tracks_roi, 100 * num_tracks_roi / num_tracks);
fprintf('  Continuous tracks: %d\n', num_continuous);
fprintf('  Discontinuous tracks: %d (rejected)\n', num_discontinuous);
fprintf('  Short continuous tracks (< %d frames): %d (rejected)\n', min_track_length_filter, num_short);
fprintf('  Long continuous tracks: %d\n', num_long + num_stationary);
fprintf('  Stationary tracks (< %d px displacement): %d (rejected)\n', min_total_displacement, num_stationary);
fprintf('  Moving tracks (>= %d px displacement): %d (KEPT)\n', min_total_displacement, num_moving);

fprintf('\nFiltering criteria:\n');
fprintf('  Max displacement/frame: %d pixels\n', max_displacement_per_frame);
fprintf('  Min track length: %d frames\n', min_track_length_filter);
fprintf('  Min total displacement: %d pixels\n', min_total_displacement);

fprintf('\nFiles saved:\n');
fprintf('  Results: %s\n', output_file);
fprintf('=============================\n');

% ========== VISUALIZE RESULTS ==========
fprintf('\nGenerating visualization...\n');

figure('Name', 'C05 - Track Filtering', 'Position', [100, 100, 1200, 800]);

% Panel 1: All tracks with short/long highlighted
subplot(1, 2, 1);
imshow(first_img);
hold on;

% Draw polygon
plot([roi_x; roi_x(1)], [roi_y; roi_y(1)], 'y-', 'LineWidth', 2);

% Show short tracks in gray
for i = 1:num_short
    track = tracks_roi_short{i};
    plot(track(:, 2), track(:, 3), '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1);
    plot(track(1, 2), track(1, 3), 'ko', 'MarkerSize', 3);
end

% Show long tracks in color
colors = jet(num_long);
for i = 1:num_long
    track = tracks_roi_long{i};
    if num_long <= 100
        plot(track(:, 2), track(:, 3), '-', 'Color', colors(i, :), 'LineWidth', 2);
    else
        plot(track(:, 2), track(:, 3), 'r-', 'LineWidth', 1);
    end
    plot(track(1, 2), track(1, 3), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
end

title(sprintf('Short (gray, N=%d) + Long (color, N=%d)', num_short, num_long));

% Panel 2: Long tracks only
subplot(1, 2, 2);
imshow(first_img);
hold on;

% Draw polygon
plot([roi_x; roi_x(1)], [roi_y; roi_y(1)], 'y-', 'LineWidth', 2);

% Show long tracks
for i = 1:num_long
    track = tracks_roi_long{i};
    if num_long <= 100
        plot(track(:, 2), track(:, 3), '-', 'Color', colors(i, :), 'LineWidth', 2);
        plot(track(1, 2), track(1, 3), 'o', 'Color', colors(i, :), ...
            'MarkerSize', 8, 'MarkerFaceColor', colors(i, :));
        plot(track(end, 2), track(end, 3), 's', 'Color', colors(i, :), ...
            'MarkerSize', 8, 'MarkerFaceColor', colors(i, :));
    else
        plot(track(:, 2), track(:, 3), 'r-', 'LineWidth', 1);
        plot(track(1, 2), track(1, 3), 'go', 'MarkerSize', 4, 'MarkerFaceColor', 'g');
        plot(track(end, 2), track(end, 3), 'rs', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    end
end

title(sprintf('Long Tracks Only (N=%d)', num_long));

fprintf('Visualization complete!\n');

fprintf('\nC05 complete!\n');
fprintf('- Polygon used: %s\n', polygon_name);
fprintf('- Discontinuous tracks rejected (frame gaps or large jumps)\n');
fprintf('- Short tracks shown in gray (rejected)\n');
fprintf('- Long continuous tracks shown in color (kept for analysis)\n');
fprintf('- Multiple polygons can be saved and selected on next run\n');
fprintf('- Next: Run C06 for Delaunay triangulation analysis\n');
