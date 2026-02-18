% C16_pca_rotation_analysis.m
% Loads filtered tracks data from C05
% Purpose: PCA rotation analysis for mantle and funnel regions
% Author: Claude
% Date: 2025-10-15

clear; close all; clc;

% ========== SELECT DATASET ==========
% Change this number to select different dataset:
%   1 = first dataset, 2 = second dataset,
DATASET_INDEX = 5;

% ========== LOAD SHARED PARAMETERS ==========
p = C00_parameters(DATASET_INDEX);

% ========== CONFIGURATION ==========
tracking_dir = p.output_dir_04;
image_dir = p.output_dir_03;
output_dir = p.output_dir_16;
image_pattern = sprintf('*.%s', p.save_format);

% Create output directory
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% ========== AUTOMATICALLY LOAD MANTLE AND FUNNEL ==========
fprintf('========== LOADING MANTLE AND FUNNEL ==========\n');

% Check for saved filtered tracks
saved_filtered_tracks = dir(fullfile(tracking_dir, 'filtered_tracks_*.mat'));

if isempty(saved_filtered_tracks)
    error('No filtered tracks found in %s\nPlease run C05 first.', tracking_dir);
end

% Find mantle and funnel files
% Pattern 1: full words "mantle" or "funnel"
% Pattern 2: abbreviated "_m_" or "_f_"
mantle_file = [];
funnel_file = [];

for i = 1:length(saved_filtered_tracks)
    fname = saved_filtered_tracks(i).name;
    fname_lower = lower(fname);

    % Check for mantle (full word or _m_)
    if contains(fname_lower, 'mantle') || contains(fname_lower, '_m_')
        mantle_file = fullfile(tracking_dir, fname);
        fprintf('Found mantle: %s\n', fname);
    % Check for funnel (full word or _f_)
    elseif contains(fname_lower, 'funnel') || contains(fname_lower, '_f_')
        funnel_file = fullfile(tracking_dir, fname);
        fprintf('Found funnel: %s\n', fname);
    end
end

if isempty(mantle_file) || isempty(funnel_file)
    fprintf('\nERROR: Could not find both mantle and funnel files.\n');
    fprintf('Available files:\n');
    for i = 1:length(saved_filtered_tracks)
        fprintf('  [%d] %s\n', i, saved_filtered_tracks(i).name);
    end
    error('Could not find both mantle and funnel files. Found: %d files', ...
        ~isempty(mantle_file) + ~isempty(funnel_file));
end

% Load mantle data
fprintf('\nLoading mantle tracks...\n');
mantle_data = load(mantle_file, 'tracks_roi_long', 'roi_x', 'roi_y', 'polygon_name');
mantle_tracks = mantle_data.tracks_roi_long;
mantle_roi_x = mantle_data.roi_x;
mantle_roi_y = mantle_data.roi_y;
if isfield(mantle_data, 'polygon_name')
    mantle_name = mantle_data.polygon_name;
else
    mantle_name = 'mantle';
end
fprintf('Loaded %d mantle tracks\n', length(mantle_tracks));

% Load funnel data
fprintf('Loading funnel tracks...\n');
funnel_data = load(funnel_file, 'tracks_roi_long', 'roi_x', 'roi_y', 'polygon_name');
funnel_tracks = funnel_data.tracks_roi_long;
funnel_roi_x = funnel_data.roi_x;
funnel_roi_y = funnel_data.roi_y;
if isfield(funnel_data, 'polygon_name')
    funnel_name = funnel_data.polygon_name;
else
    funnel_name = 'funnel';
end
fprintf('Loaded %d funnel tracks\n', length(funnel_tracks));

% ========== PROCESS MANTLE TRACKS ==========
fprintf('\n========== PROCESSING MANTLE ==========\n');
[mantle_tracks, mantle_start_pos, mantle_num_full, mantle_num_partial] = ...
    process_tracks(mantle_tracks, p.require_full_duration_tracks);

% ========== PROCESS FUNNEL TRACKS ==========
fprintf('\n========== PROCESSING FUNNEL ==========\n');
[funnel_tracks, funnel_start_pos, funnel_num_full, funnel_num_partial] = ...
    process_tracks(funnel_tracks, p.require_full_duration_tracks);

% Determine number of frames from both regions
mantle_lengths = cellfun(@(x) size(x, 1), mantle_tracks);
funnel_lengths = cellfun(@(x) size(x, 1), funnel_tracks);
num_frames = max([max(mantle_lengths), max(funnel_lengths)]);

% Limit frames if specified in parameters
if p.num_frames_to_analyze > 0 && p.num_frames_to_analyze < num_frames
    num_frames = p.num_frames_to_analyze;
    fprintf('\nLimiting to first %d frames (from parameter)\n', num_frames);
else
    fprintf('\nTotal frames: %d\n', num_frames);
end

% Load first frame
image_files = dir(fullfile(image_dir, image_pattern));
first_img = imread(fullfile(image_dir, image_files(1).name));

% ========== COMPUTE PCA FOR MANTLE ==========
fprintf('\n========== COMPUTING PCA FOR MANTLE ==========\n');
[mantle_pca_coeff1, mantle_pca_coeff2, mantle_pca_center, mantle_pca_explained, mantle_pca_eigenvalues] = ...
    compute_pca_all_frames(mantle_tracks, num_frames);

fprintf('Mantle PCA Statistics:\n');
fprintf('  Mean explained variance (PC1): %.2f%%\n', mean(mantle_pca_explained(:, 1)));
fprintf('  Mean explained variance (PC2): %.2f%%\n', mean(mantle_pca_explained(:, 2)));
fprintf('  Mean eigenvalue ratio (PC1/PC2): %.2f\n', ...
    mean(mantle_pca_eigenvalues(:, 1) ./ mantle_pca_eigenvalues(:, 2)));

% ========== COMPUTE PCA FOR FUNNEL ==========
fprintf('\n========== COMPUTING PCA FOR FUNNEL ==========\n');
[funnel_pca_coeff1, funnel_pca_coeff2, funnel_pca_center, funnel_pca_explained, funnel_pca_eigenvalues] = ...
    compute_pca_all_frames(funnel_tracks, num_frames);

fprintf('Funnel PCA Statistics:\n');
fprintf('  Mean explained variance (PC1): %.2f%%\n', mean(funnel_pca_explained(:, 1)));
fprintf('  Mean explained variance (PC2): %.2f%%\n', mean(funnel_pca_explained(:, 2)));
fprintf('  Mean eigenvalue ratio (PC1/PC2): %.2f\n', ...
    mean(funnel_pca_eigenvalues(:, 1) ./ funnel_pca_eigenvalues(:, 2)));

% ========== COMPUTE ROTATED LAB COORDINATE PROJECTIONS ==========
fprintf('\n========== COMPUTING ROTATED LAB COORDINATE PROJECTIONS ==========\n');

% Mantle analysis
fprintf('Analyzing mantle...\n');
[mantle_avg_dist_x, mantle_avg_dist_y, mantle_rotation_angle] = ...
    compute_rotated_lab_distances(mantle_tracks, mantle_pca_coeff1, mantle_pca_center, num_frames);

% Funnel analysis
fprintf('Analyzing funnel...\n');
[funnel_avg_dist_x, funnel_avg_dist_y, funnel_rotation_angle] = ...
    compute_rotated_lab_distances(funnel_tracks, funnel_pca_coeff1, funnel_pca_center, num_frames);

fprintf('Rotated lab coordinate analysis complete!\n');

% ========== VISUALIZE PCA RESULTS ==========
fprintf('\nGenerating visualizations...\n');

% Convert frames to time
time_ms = (1:num_frames) / p.frame_rate * 1000;

% Load raw images
raw_image_dir = p.data_dir;
raw_image_files = dir(fullfile(raw_image_dir, p.image_pattern));

% Select three time points: start, middle, end
frames_to_show = [1, round(num_frames/2), num_frames];

% Figure 1: Combined Mantle and Funnel PCA rotation over time (overlaid)
figure('Name', 'C16 - Mantle & Funnel PCA Rotation', 'Position', [50, 100, 1800, 600]);
visualize_pca_rotation_combined(raw_image_files, raw_image_dir, frames_to_show, time_ms, ...
    mantle_tracks, mantle_roi_x, mantle_roi_y, mantle_pca_coeff1, mantle_pca_coeff2, ...
    mantle_pca_center, mantle_pca_eigenvalues, mantle_rotation_angle, ...
    funnel_tracks, funnel_roi_x, funnel_roi_y, funnel_pca_coeff1, funnel_pca_coeff2, ...
    funnel_pca_center, funnel_pca_eigenvalues, funnel_rotation_angle);
saveas(gcf, fullfile(output_dir, 'C16_combined_pca_rotation.png'));

% Figure 2: Rotated lab coordinate distances
figure('Name', 'C16 - Rotated Lab Distances', 'Position', [50, 100, 1600, 800]);

% Mantle distances
subplot(2, 3, 1);
plot(time_ms, mantle_avg_dist_x, 'r-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Average Distance (pixels)');
title('Mantle - Along PC1 (Tissue X)');
grid on;

subplot(2, 3, 2);
plot(time_ms, mantle_avg_dist_y, 'b-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Average Distance (pixels)');
title('Mantle - Perpendicular to PC1 (Tissue Y)');
grid on;

subplot(2, 3, 3);
plot(time_ms, mantle_rotation_angle, 'k-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Rotation Angle (degrees)');
title('Mantle - Rotation from Initial');
grid on;

% Funnel distances
subplot(2, 3, 4);
plot(time_ms, funnel_avg_dist_x, 'r-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Average Distance (pixels)');
title('Funnel - Along PC1 (Tissue X)');
grid on;

subplot(2, 3, 5);
plot(time_ms, funnel_avg_dist_y, 'b-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Average Distance (pixels)');
title('Funnel - Perpendicular to PC1 (Tissue Y)');
grid on;

subplot(2, 3, 6);
plot(time_ms, funnel_rotation_angle, 'k-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Rotation Angle (degrees)');
title('Funnel - Rotation from Initial');
grid on;

saveas(gcf, fullfile(output_dir, 'C16_rotated_lab_distances.png'));

fprintf('Figures saved!\n');

% ========== SAVE RESULTS ==========
fprintf('\nSaving results...\n');

output_file = fullfile(output_dir, 'C16_pca_results_both.mat');

save(output_file, ...
    'mantle_pca_coeff1', 'mantle_pca_coeff2', 'mantle_pca_center', ...
    'mantle_pca_explained', 'mantle_pca_eigenvalues', 'mantle_tracks', ...
    'mantle_roi_x', 'mantle_roi_y', 'mantle_name', ...
    'mantle_avg_dist_x', 'mantle_avg_dist_y', 'mantle_rotation_angle', ...
    'funnel_pca_coeff1', 'funnel_pca_coeff2', 'funnel_pca_center', ...
    'funnel_pca_explained', 'funnel_pca_eigenvalues', 'funnel_tracks', ...
    'funnel_roi_x', 'funnel_roi_y', 'funnel_name', ...
    'funnel_avg_dist_x', 'funnel_avg_dist_y', 'funnel_rotation_angle', ...
    'num_frames');

fprintf('Saved to: %s\n', output_file);

% ========== ASK USER FOR VIDEO GENERATION ==========
fprintf('\n');
user_response = input('Generate video with PCA overlay for both regions? (This will take time) (y/n): ', 's');

if strcmpi(user_response, 'y')
    % Generate combined video showing both regions
    fprintf('\n========== GENERATING COMBINED VIDEO ==========\n');
    generate_combined_pca_video(output_dir, raw_image_dir, raw_image_files, num_frames, ...
        mantle_tracks, mantle_roi_x, mantle_roi_y, mantle_pca_coeff1, mantle_pca_coeff2, ...
        mantle_pca_center, mantle_pca_eigenvalues, mantle_rotation_angle, ...
        funnel_tracks, funnel_roi_x, funnel_roi_y, funnel_pca_coeff1, funnel_pca_coeff2, ...
        funnel_pca_center, funnel_pca_eigenvalues, funnel_rotation_angle, ...
        p.frame_rate);
else
    fprintf('Skipping video generation.\n');
end

% ========== PRINT SUMMARY ==========
fprintf('\n========== SUMMARY ==========\n');
fprintf('Mantle:\n');
fprintf('  Tracks loaded: %d\n', mantle_num_full + mantle_num_partial);
fprintf('  Full duration: %d (used)\n', mantle_num_full);
fprintf('  Partial duration: %d (rejected)\n', mantle_num_partial);
fprintf('  Mean elongation: %.2f\n', mean(mantle_pca_eigenvalues(:, 1) ./ mantle_pca_eigenvalues(:, 2)));

fprintf('\nFunnel:\n');
fprintf('  Tracks loaded: %d\n', funnel_num_full + funnel_num_partial);
fprintf('  Full duration: %d (used)\n', funnel_num_full);
fprintf('  Partial duration: %d (rejected)\n', funnel_num_partial);
fprintf('  Mean elongation: %.2f\n', mean(funnel_pca_eigenvalues(:, 1) ./ funnel_pca_eigenvalues(:, 2)));

fprintf('\nGeneral:\n');
fprintf('  Total frames: %d\n', num_frames);
fprintf('  Frame rate: %.2f fps\n', p.frame_rate);
fprintf('  Duration: %.2f seconds\n', num_frames / p.frame_rate);
fprintf('=============================\n');

fprintf('\nC16 complete!\n');
fprintf('- PCA computed for both mantle and funnel\n');
fprintf('- Results saved to: %s\n', output_file);

%% ========== HELPER FUNCTIONS ==========

function [tracks_out, start_pos, num_full, num_partial] = process_tracks(tracks_in, require_full_duration)
    % Process tracks: filter by duration and remove duplicates

    if require_full_duration
        fprintf('Filtering tracks to keep only those with full duration...\n');

        % Determine maximum duration
        track_lengths = cellfun(@(x) size(x, 1), tracks_in);
        max_duration = max(track_lengths);
        fprintf('Maximum track duration: %d frames\n', max_duration);

        % Keep only tracks that cover full duration
        full_duration_mask = track_lengths == max_duration;
        tracks_full = tracks_in(full_duration_mask);

        num_full = sum(full_duration_mask);
        num_partial = sum(~full_duration_mask);

        fprintf('Full duration tracks: %d (KEPT)\n', num_full);
        fprintf('Partial duration tracks: %d (rejected)\n', num_partial);

        tracks_in = tracks_full;
    else
        fprintf('Using all filtered tracks (full duration filter disabled)...\n');
        num_full = length(tracks_in);
        num_partial = 0;
    end

    % Get starting positions
    num_tracks = length(tracks_in);
    start_positions = zeros(num_tracks, 2);
    for i = 1:num_tracks
        track = tracks_in{i};
        start_positions(i, :) = track(1, 2:3);
    end

    % Remove duplicate starting positions
    [unique_positions, unique_idx, ~] = unique(start_positions, 'rows', 'stable');
    num_duplicates = num_tracks - size(unique_positions, 1);

    if num_duplicates > 0
        fprintf('WARNING: Found %d duplicate starting positions - removing duplicates\n', num_duplicates);
        tracks_in = tracks_in(unique_idx);
        start_positions = unique_positions;
    end

    tracks_out = tracks_in;
    start_pos = start_positions;
    fprintf('Final track count: %d\n', length(tracks_out));
end

function [coeff1, coeff2, center, explained, eigenvalues] = compute_pca_all_frames(tracks, num_frames)
    % Compute PCA for all frames

    num_tracks = length(tracks);

    % Preallocate arrays for PCA results
    coeff1 = zeros(num_frames, 2);  % First principal component (major axis)
    coeff2 = zeros(num_frames, 2);  % Second principal component (minor axis)
    center = zeros(num_frames, 2);  % Center (mean position)
    explained = zeros(num_frames, 2);  % Explained variance
    eigenvalues = zeros(num_frames, 2);  % Eigenvalues

    % Compute PCA for each frame
    for frame = 1:num_frames
        % Get positions at this frame
        positions = [];
        for i = 1:num_tracks
            track = tracks{i};
            if frame <= size(track, 1)
                positions = [positions; track(frame, 2:3)];
            end
        end

        % Need at least 2 points for PCA
        if size(positions, 1) >= 2
            % Center the data
            mean_pos = mean(positions, 1);
            centered_data = positions - mean_pos;

            % Compute covariance matrix
            cov_matrix = cov(centered_data);

            % Compute eigenvalues and eigenvectors
            [eigenvectors, eigenval_matrix] = eig(cov_matrix);
            eigenval_diag = diag(eigenval_matrix);

            % Sort by eigenvalues (largest first)
            [eigenval_sorted, idx] = sort(eigenval_diag, 'descend');
            eigenvec_sorted = eigenvectors(:, idx);

            % Store results
            center(frame, :) = mean_pos;
            coeff1(frame, :) = eigenvec_sorted(:, 1)';  % Major axis
            coeff2(frame, :) = eigenvec_sorted(:, 2)';  % Minor axis
            eigenvalues(frame, :) = eigenval_sorted';

            % Explained variance
            total_variance = sum(eigenval_sorted);
            if total_variance > 0
                explained(frame, :) = 100 * eigenval_sorted' / total_variance;
            end
        end

        if mod(frame, 50) == 0 || frame == num_frames
            fprintf('  Processed frame %d/%d\n', frame, num_frames);
        end
    end

    fprintf('PCA computation complete!\n');
end

function visualize_pca_rotation(image_files, image_dir, frames_to_show, time_ms, tracks, roi_x, roi_y, pca_coeff1, pca_coeff2, pca_center, pca_eigenvalues, region_name)
    % Visualize PCA rotation at three time points

    for i = 1:3
        frame = frames_to_show(i);

        % Load image for this frame
        if frame <= length(image_files)
            img = imread(fullfile(image_dir, image_files(frame).name));
        else
            img = imread(fullfile(image_dir, image_files(end).name));
        end

        subplot(1, 3, i);
        imshow(img);
        hold on;

        % Draw polygon
        plot([roi_x; roi_x(1)], [roi_y; roi_y(1)], 'y-', 'LineWidth', 2);

        % Get positions at this frame
        positions = [];
        for j = 1:length(tracks)
            track = tracks{j};
            if frame <= size(track, 1)
                positions = [positions; track(frame, 2:3)];
            end
        end

        % Draw points
        if ~isempty(positions)
            plot(positions(:, 1), positions(:, 2), 'ro', ...
                'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
        end

        % Draw PCA axes
        center_pos = pca_center(frame, :);
        pc1 = pca_coeff1(frame, :);
        pc2 = pca_coeff2(frame, :);
        scale = 3 * sqrt(pca_eigenvalues(frame, 1));

        % Major axis (PC1) - red
        quiver(center_pos(1), center_pos(2), pc1(1)*scale, pc1(2)*scale, 0, ...
            'r-', 'LineWidth', 3, 'MaxHeadSize', 2);
        quiver(center_pos(1), center_pos(2), -pc1(1)*scale, -pc1(2)*scale, 0, ...
            'r-', 'LineWidth', 3, 'MaxHeadSize', 2);

        % Minor axis (PC2) - cyan
        scale2 = scale * sqrt(pca_eigenvalues(frame, 2) / pca_eigenvalues(frame, 1));
        quiver(center_pos(1), center_pos(2), pc2(1)*scale2, pc2(2)*scale2, 0, ...
            'c-', 'LineWidth', 3, 'MaxHeadSize', 2);
        quiver(center_pos(1), center_pos(2), -pc2(1)*scale2, -pc2(2)*scale2, 0, ...
            'c-', 'LineWidth', 3, 'MaxHeadSize', 2);

        % Draw center
        plot(center_pos(1), center_pos(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y', 'LineWidth', 2);

        % Calculate angle and elongation
        angle = atan2d(pc1(2), pc1(1));
        elongation = pca_eigenvalues(frame, 1) / pca_eigenvalues(frame, 2);

        title(sprintf('%s - Frame %d (%.1f ms)\nAngle: %.1f°, Elongation: %.2f', ...
            region_name, frame, time_ms(frame), angle, elongation), 'FontSize', 10);
    end
end

function visualize_pca_rotation_combined(image_files, image_dir, frames_to_show, time_ms, ...
    mantle_tracks, mantle_roi_x, mantle_roi_y, mantle_pca_coeff1, mantle_pca_coeff2, ...
    mantle_pca_center, mantle_pca_eigenvalues, mantle_rotation_angle, ...
    funnel_tracks, funnel_roi_x, funnel_roi_y, funnel_pca_coeff1, funnel_pca_coeff2, ...
    funnel_pca_center, funnel_pca_eigenvalues, funnel_rotation_angle)
    % Visualize both mantle and funnel PCA rotation overlaid on same images
    % Also shows rotated lab coordinate axes

    for i = 1:3
        frame = frames_to_show(i);

        % Load image for this frame
        if frame <= length(image_files)
            img = imread(fullfile(image_dir, image_files(frame).name));
        else
            img = imread(fullfile(image_dir, image_files(end).name));
        end

        subplot(1, 3, i);
        imshow(img);
        hold on;

        % Draw mantle polygon
        plot([mantle_roi_x; mantle_roi_x(1)], [mantle_roi_y; mantle_roi_y(1)], ...
            'r-', 'LineWidth', 2);

        % Draw funnel polygon
        plot([funnel_roi_x; funnel_roi_x(1)], [funnel_roi_y; funnel_roi_y(1)], ...
            'b-', 'LineWidth', 2);

        % Get mantle positions at this frame
        mantle_positions = [];
        for j = 1:length(mantle_tracks)
            track = mantle_tracks{j};
            if frame <= size(track, 1)
                mantle_positions = [mantle_positions; track(frame, 2:3)];
            end
        end

        % Get funnel positions at this frame
        funnel_positions = [];
        for j = 1:length(funnel_tracks)
            track = funnel_tracks{j};
            if frame <= size(track, 1)
                funnel_positions = [funnel_positions; track(frame, 2:3)];
            end
        end

        % Draw mantle points (red)
        if ~isempty(mantle_positions)
            plot(mantle_positions(:, 1), mantle_positions(:, 2), 'ro', ...
                'MarkerSize', 5, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'r', 'LineWidth', 1);
        end

        % Draw funnel points (blue)
        if ~isempty(funnel_positions)
            plot(funnel_positions(:, 1), funnel_positions(:, 2), 'bo', ...
                'MarkerSize', 5, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b', 'LineWidth', 1);
        end

        % Draw mantle PCA axes
        mantle_center = mantle_pca_center(frame, :);
        mantle_pc1 = mantle_pca_coeff1(frame, :);
        mantle_pc2 = mantle_pca_coeff2(frame, :);
        mantle_scale = 3 * sqrt(mantle_pca_eigenvalues(frame, 1));

        % Mantle PC1 - red
        quiver(mantle_center(1), mantle_center(2), mantle_pc1(1)*mantle_scale, mantle_pc1(2)*mantle_scale, 0, ...
            'r-', 'LineWidth', 3, 'MaxHeadSize', 2);
        quiver(mantle_center(1), mantle_center(2), -mantle_pc1(1)*mantle_scale, -mantle_pc1(2)*mantle_scale, 0, ...
            'r-', 'LineWidth', 3, 'MaxHeadSize', 2);
        % PC1 label
        text_pos = mantle_center + 1.2 * mantle_scale * mantle_pc1;
        text(text_pos(1), text_pos(2), 'PC1', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');

        % Mantle PC2 - magenta
        mantle_scale2 = mantle_scale * sqrt(mantle_pca_eigenvalues(frame, 2) / mantle_pca_eigenvalues(frame, 1));
        quiver(mantle_center(1), mantle_center(2), mantle_pc2(1)*mantle_scale2, mantle_pc2(2)*mantle_scale2, 0, ...
            'm-', 'LineWidth', 2, 'MaxHeadSize', 2);
        quiver(mantle_center(1), mantle_center(2), -mantle_pc2(1)*mantle_scale2, -mantle_pc2(2)*mantle_scale2, 0, ...
            'm-', 'LineWidth', 2, 'MaxHeadSize', 2);

        % Mantle center
        plot(mantle_center(1), mantle_center(2), 'ko', 'MarkerSize', 10, ...
            'MarkerFaceColor', 'r', 'LineWidth', 2);

        % Draw funnel PCA axes
        funnel_center = funnel_pca_center(frame, :);
        funnel_pc1 = funnel_pca_coeff1(frame, :);
        funnel_pc2 = funnel_pca_coeff2(frame, :);
        funnel_scale = 3 * sqrt(funnel_pca_eigenvalues(frame, 1));

        % Funnel PC1 - blue
        quiver(funnel_center(1), funnel_center(2), funnel_pc1(1)*funnel_scale, funnel_pc1(2)*funnel_scale, 0, ...
            'b-', 'LineWidth', 3, 'MaxHeadSize', 2);
        quiver(funnel_center(1), funnel_center(2), -funnel_pc1(1)*funnel_scale, -funnel_pc1(2)*funnel_scale, 0, ...
            'b-', 'LineWidth', 3, 'MaxHeadSize', 2);
        % PC1 label
        text_pos = funnel_center + 1.2 * funnel_scale * funnel_pc1;
        text(text_pos(1), text_pos(2), 'PC1', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');

        % Funnel PC2 - cyan
        funnel_scale2 = funnel_scale * sqrt(funnel_pca_eigenvalues(frame, 2) / funnel_pca_eigenvalues(frame, 1));
        quiver(funnel_center(1), funnel_center(2), funnel_pc2(1)*funnel_scale2, funnel_pc2(2)*funnel_scale2, 0, ...
            'c-', 'LineWidth', 2, 'MaxHeadSize', 2);
        quiver(funnel_center(1), funnel_center(2), -funnel_pc2(1)*funnel_scale2, -funnel_pc2(2)*funnel_scale2, 0, ...
            'c-', 'LineWidth', 2, 'MaxHeadSize', 2);

        % Funnel center
        plot(funnel_center(1), funnel_center(2), 'ko', 'MarkerSize', 10, ...
            'MarkerFaceColor', 'b', 'LineWidth', 2);

        % ===== DRAW ROTATED LAB COORDINATES =====
        % Lab starts aligned with image axes at frame 1
        % Then rotates by the relative angle change of PCA
        lab_scale = 60;  % Fixed length for lab axes

        % Get initial and current PCA angles
        mantle_initial_angle = atan2(mantle_pca_coeff1(1, 2), mantle_pca_coeff1(1, 1));
        mantle_current_angle = atan2(mantle_pc1(2), mantle_pc1(1));
        funnel_initial_angle = atan2(funnel_pca_coeff1(1, 2), funnel_pca_coeff1(1, 1));
        funnel_current_angle = atan2(funnel_pc1(2), funnel_pc1(1));

        % Lab angle = initial (0 = horizontal) + relative rotation
        mantle_lab_angle = mantle_current_angle - mantle_initial_angle;
        funnel_lab_angle = funnel_current_angle - funnel_initial_angle;

        % Mantle lab axes (dashed yellow)
        mantle_lab_x = [cos(mantle_lab_angle), sin(mantle_lab_angle)];
        mantle_lab_y = [-sin(mantle_lab_angle), cos(mantle_lab_angle)];

        quiver(mantle_center(1), mantle_center(2), mantle_lab_x(1)*lab_scale, mantle_lab_x(2)*lab_scale, 0, ...
            'y--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(mantle_center(1), mantle_center(2), -mantle_lab_x(1)*lab_scale, -mantle_lab_x(2)*lab_scale, 0, ...
            'y--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(mantle_center(1), mantle_center(2), mantle_lab_y(1)*lab_scale, mantle_lab_y(2)*lab_scale, 0, ...
            'y--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(mantle_center(1), mantle_center(2), -mantle_lab_y(1)*lab_scale, -mantle_lab_y(2)*lab_scale, 0, ...
            'y--', 'LineWidth', 2, 'MaxHeadSize', 1.5);

        % Funnel lab axes (dashed green)
        funnel_lab_x = [cos(funnel_lab_angle), sin(funnel_lab_angle)];
        funnel_lab_y = [-sin(funnel_lab_angle), cos(funnel_lab_angle)];

        quiver(funnel_center(1), funnel_center(2), funnel_lab_x(1)*lab_scale, funnel_lab_x(2)*lab_scale, 0, ...
            'g--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(funnel_center(1), funnel_center(2), -funnel_lab_x(1)*lab_scale, -funnel_lab_x(2)*lab_scale, 0, ...
            'g--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(funnel_center(1), funnel_center(2), funnel_lab_y(1)*lab_scale, funnel_lab_y(2)*lab_scale, 0, ...
            'g--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(funnel_center(1), funnel_center(2), -funnel_lab_y(1)*lab_scale, -funnel_lab_y(2)*lab_scale, 0, ...
            'g--', 'LineWidth', 2, 'MaxHeadSize', 1.5);

        % ===== DRAW PROJECTIONS ONTO VERTICAL AXES =====
        % Project mantle points onto mantle lab Y-axis
        if ~isempty(mantle_positions)
            for j = 1:size(mantle_positions, 1)
                pos = mantle_positions(j, :);
                % Vector from center to point
                vec = pos - mantle_center;
                % Projection onto lab Y-axis
                proj_length = dot(vec, mantle_lab_y);
                proj_point = mantle_center + proj_length * mantle_lab_y;
                % Draw projection line
                plot([pos(1), proj_point(1)], [pos(2), proj_point(2)], ...
                    'r:', 'LineWidth', 0.5);
                % Draw projected point
                plot(proj_point(1), proj_point(2), 'r^', ...
                    'MarkerSize', 4, 'MarkerFaceColor', 'yellow');
            end
        end

        % Project funnel points onto funnel lab Y-axis
        if ~isempty(funnel_positions)
            for j = 1:size(funnel_positions, 1)
                pos = funnel_positions(j, :);
                % Vector from center to point
                vec = pos - funnel_center;
                % Projection onto lab Y-axis
                proj_length = dot(vec, funnel_lab_y);
                proj_point = funnel_center + proj_length * funnel_lab_y;
                % Draw projection line
                plot([pos(1), proj_point(1)], [pos(2), proj_point(2)], ...
                    'b:', 'LineWidth', 0.5);
                % Draw projected point
                plot(proj_point(1), proj_point(2), 'b^', ...
                    'MarkerSize', 4, 'MarkerFaceColor', 'cyan');
            end
        end

        % Calculate angles and elongations
        mantle_angle = atan2d(mantle_pc1(2), mantle_pc1(1));
        mantle_elong = mantle_pca_eigenvalues(frame, 1) / mantle_pca_eigenvalues(frame, 2);
        funnel_angle = atan2d(funnel_pc1(2), funnel_pc1(1));
        funnel_elong = funnel_pca_eigenvalues(frame, 1) / funnel_pca_eigenvalues(frame, 2);

        title(sprintf('Frame %d (%.1f ms)\nMantle: Δ%.1f°, Elong=%.2f | Funnel: Δ%.1f°, Elong=%.2f', ...
            frame, time_ms(frame), mantle_rotation_angle(frame), mantle_elong, ...
            funnel_rotation_angle(frame), funnel_elong), 'FontSize', 9);
    end
end

function generate_combined_pca_video(output_dir, image_dir, image_files, num_frames, ...
    mantle_tracks, mantle_roi_x, mantle_roi_y, mantle_pca_coeff1, mantle_pca_coeff2, ...
    mantle_pca_center, mantle_pca_eigenvalues, mantle_rotation_angle, ...
    funnel_tracks, funnel_roi_x, funnel_roi_y, funnel_pca_coeff1, funnel_pca_coeff2, ...
    funnel_pca_center, funnel_pca_eigenvalues, funnel_rotation_angle, frame_rate)
    % Generate video with both mantle and funnel PCA overlays (exactly like Figure 1)

    % Create output subfolder
    frames_dir = fullfile(output_dir, 'C16_pca_frames_combined');
    if ~exist(frames_dir, 'dir')
        mkdir(frames_dir);
    end

    % Downsample frames for visualization
    if num_frames > 100
        frame_step = 5;
    else
        frame_step = 1;
    end

    frames_to_render = 1:frame_step:num_frames;
    fprintf('Rendering %d frames (every %d frame)...\n', length(frames_to_render), frame_step);

    % Initialize video writer
    video_filename = fullfile(output_dir, 'C16_pca_video_combined.mp4');
    video_writer = [];

    % Lab scale
    lab_scale = 60;

    for idx = 1:length(frames_to_render)
        frame = frames_to_render(idx);

        % Load raw image
        if frame <= length(image_files)
            img = imread(fullfile(image_dir, image_files(frame).name));
        else
            img = imread(fullfile(image_dir, image_files(end).name));
        end

        % Get image size
        [img_height, img_width, ~] = size(img);

        % Create figure
        fig = figure('Visible', 'off', 'Units', 'pixels', ...
            'Position', [100, 100, img_width, img_height]);
        ax = axes('Units', 'pixels', 'Position', [0, 0, img_width, img_height]);
        imshow(img, 'Parent', ax);
        hold on;

        % Draw mantle polygon
        plot([mantle_roi_x; mantle_roi_x(1)], [mantle_roi_y; mantle_roi_y(1)], ...
            'r-', 'LineWidth', 2);

        % Draw funnel polygon
        plot([funnel_roi_x; funnel_roi_x(1)], [funnel_roi_y; funnel_roi_y(1)], ...
            'b-', 'LineWidth', 2);

        % Get mantle positions at this frame
        mantle_positions = [];
        for j = 1:length(mantle_tracks)
            track = mantle_tracks{j};
            if frame <= size(track, 1)
                mantle_positions = [mantle_positions; track(frame, 2:3)];
            end
        end

        % Get funnel positions at this frame
        funnel_positions = [];
        for j = 1:length(funnel_tracks)
            track = funnel_tracks{j};
            if frame <= size(track, 1)
                funnel_positions = [funnel_positions; track(frame, 2:3)];
            end
        end

        % Draw mantle points (red)
        if ~isempty(mantle_positions)
            plot(mantle_positions(:, 1), mantle_positions(:, 2), 'ro', ...
                'MarkerSize', 5, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'r', 'LineWidth', 1);
        end

        % Draw funnel points (blue)
        if ~isempty(funnel_positions)
            plot(funnel_positions(:, 1), funnel_positions(:, 2), 'bo', ...
                'MarkerSize', 5, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b', 'LineWidth', 1);
        end

        % Draw mantle PCA axes
        mantle_center = mantle_pca_center(frame, :);
        mantle_pc1 = mantle_pca_coeff1(frame, :);
        mantle_pc2 = mantle_pca_coeff2(frame, :);
        mantle_scale = 3 * sqrt(mantle_pca_eigenvalues(frame, 1));

        % Mantle PC1 - red
        quiver(mantle_center(1), mantle_center(2), mantle_pc1(1)*mantle_scale, mantle_pc1(2)*mantle_scale, 0, ...
            'r-', 'LineWidth', 3, 'MaxHeadSize', 2);
        quiver(mantle_center(1), mantle_center(2), -mantle_pc1(1)*mantle_scale, -mantle_pc1(2)*mantle_scale, 0, ...
            'r-', 'LineWidth', 3, 'MaxHeadSize', 2);
        % PC1 label
        text_pos = mantle_center + 1.2 * mantle_scale * mantle_pc1;
        text(text_pos(1), text_pos(2), 'PC1', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');

        % Mantle PC2 - magenta
        mantle_scale2 = mantle_scale * sqrt(mantle_pca_eigenvalues(frame, 2) / mantle_pca_eigenvalues(frame, 1));
        quiver(mantle_center(1), mantle_center(2), mantle_pc2(1)*mantle_scale2, mantle_pc2(2)*mantle_scale2, 0, ...
            'm-', 'LineWidth', 2, 'MaxHeadSize', 2);
        quiver(mantle_center(1), mantle_center(2), -mantle_pc2(1)*mantle_scale2, -mantle_pc2(2)*mantle_scale2, 0, ...
            'm-', 'LineWidth', 2, 'MaxHeadSize', 2);

        % Mantle center
        plot(mantle_center(1), mantle_center(2), 'ko', 'MarkerSize', 10, ...
            'MarkerFaceColor', 'r', 'LineWidth', 2);

        % Draw funnel PCA axes
        funnel_center = funnel_pca_center(frame, :);
        funnel_pc1 = funnel_pca_coeff1(frame, :);
        funnel_pc2 = funnel_pca_coeff2(frame, :);
        funnel_scale = 3 * sqrt(funnel_pca_eigenvalues(frame, 1));

        % Funnel PC1 - blue
        quiver(funnel_center(1), funnel_center(2), funnel_pc1(1)*funnel_scale, funnel_pc1(2)*funnel_scale, 0, ...
            'b-', 'LineWidth', 3, 'MaxHeadSize', 2);
        quiver(funnel_center(1), funnel_center(2), -funnel_pc1(1)*funnel_scale, -funnel_pc1(2)*funnel_scale, 0, ...
            'b-', 'LineWidth', 3, 'MaxHeadSize', 2);
        % PC1 label
        text_pos = funnel_center + 1.2 * funnel_scale * funnel_pc1;
        text(text_pos(1), text_pos(2), 'PC1', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');

        % Funnel PC2 - cyan
        funnel_scale2 = funnel_scale * sqrt(funnel_pca_eigenvalues(frame, 2) / funnel_pca_eigenvalues(frame, 1));
        quiver(funnel_center(1), funnel_center(2), funnel_pc2(1)*funnel_scale2, funnel_pc2(2)*funnel_scale2, 0, ...
            'c-', 'LineWidth', 2, 'MaxHeadSize', 2);
        quiver(funnel_center(1), funnel_center(2), -funnel_pc2(1)*funnel_scale2, -funnel_pc2(2)*funnel_scale2, 0, ...
            'c-', 'LineWidth', 2, 'MaxHeadSize', 2);

        % Funnel center
        plot(funnel_center(1), funnel_center(2), 'ko', 'MarkerSize', 10, ...
            'MarkerFaceColor', 'b', 'LineWidth', 2);

        % ===== DRAW ROTATED LAB COORDINATES =====
        % Get initial and current PCA angles
        mantle_initial_angle = atan2(mantle_pca_coeff1(1, 2), mantle_pca_coeff1(1, 1));
        mantle_current_angle = atan2(mantle_pc1(2), mantle_pc1(1));
        funnel_initial_angle = atan2(funnel_pca_coeff1(1, 2), funnel_pca_coeff1(1, 1));
        funnel_current_angle = atan2(funnel_pc1(2), funnel_pc1(1));

        % Lab angle = initial (0 = horizontal) + relative rotation
        mantle_lab_angle = mantle_current_angle - mantle_initial_angle;
        funnel_lab_angle = funnel_current_angle - funnel_initial_angle;

        % Mantle lab axes (dashed yellow)
        mantle_lab_x = [cos(mantle_lab_angle), sin(mantle_lab_angle)];
        mantle_lab_y = [-sin(mantle_lab_angle), cos(mantle_lab_angle)];

        quiver(mantle_center(1), mantle_center(2), mantle_lab_x(1)*lab_scale, mantle_lab_x(2)*lab_scale, 0, ...
            'y--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(mantle_center(1), mantle_center(2), -mantle_lab_x(1)*lab_scale, -mantle_lab_x(2)*lab_scale, 0, ...
            'y--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(mantle_center(1), mantle_center(2), mantle_lab_y(1)*lab_scale, mantle_lab_y(2)*lab_scale, 0, ...
            'y--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(mantle_center(1), mantle_center(2), -mantle_lab_y(1)*lab_scale, -mantle_lab_y(2)*lab_scale, 0, ...
            'y--', 'LineWidth', 2, 'MaxHeadSize', 1.5);

        % Funnel lab axes (dashed green)
        funnel_lab_x = [cos(funnel_lab_angle), sin(funnel_lab_angle)];
        funnel_lab_y = [-sin(funnel_lab_angle), cos(funnel_lab_angle)];

        quiver(funnel_center(1), funnel_center(2), funnel_lab_x(1)*lab_scale, funnel_lab_x(2)*lab_scale, 0, ...
            'g--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(funnel_center(1), funnel_center(2), -funnel_lab_x(1)*lab_scale, -funnel_lab_x(2)*lab_scale, 0, ...
            'g--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(funnel_center(1), funnel_center(2), funnel_lab_y(1)*lab_scale, funnel_lab_y(2)*lab_scale, 0, ...
            'g--', 'LineWidth', 2, 'MaxHeadSize', 1.5);
        quiver(funnel_center(1), funnel_center(2), -funnel_lab_y(1)*lab_scale, -funnel_lab_y(2)*lab_scale, 0, ...
            'g--', 'LineWidth', 2, 'MaxHeadSize', 1.5);

        % ===== DRAW PROJECTIONS ONTO VERTICAL AXES =====
        % Project mantle points onto mantle lab Y-axis
        if ~isempty(mantle_positions)
            for j = 1:size(mantle_positions, 1)
                pos = mantle_positions(j, :);
                vec = pos - mantle_center;
                proj_length = dot(vec, mantle_lab_y);
                proj_point = mantle_center + proj_length * mantle_lab_y;
                plot([pos(1), proj_point(1)], [pos(2), proj_point(2)], 'r:', 'LineWidth', 0.5);
                plot(proj_point(1), proj_point(2), 'r^', 'MarkerSize', 4, 'MarkerFaceColor', 'yellow');
            end
        end

        % Project funnel points onto funnel lab Y-axis
        if ~isempty(funnel_positions)
            for j = 1:size(funnel_positions, 1)
                pos = funnel_positions(j, :);
                vec = pos - funnel_center;
                proj_length = dot(vec, funnel_lab_y);
                proj_point = funnel_center + proj_length * funnel_lab_y;
                plot([pos(1), proj_point(1)], [pos(2), proj_point(2)], 'b:', 'LineWidth', 0.5);
                plot(proj_point(1), proj_point(2), 'b^', 'MarkerSize', 4, 'MarkerFaceColor', 'cyan');
            end
        end

        % Calculate angles and elongations
        mantle_elong = mantle_pca_eigenvalues(frame, 1) / mantle_pca_eigenvalues(frame, 2);
        funnel_elong = funnel_pca_eigenvalues(frame, 1) / funnel_pca_eigenvalues(frame, 2);

        % Add title with frame info
        time_ms = (frame - 1) * (1000 / frame_rate);
        title(sprintf('Frame %d (%.1f ms) | Mantle: Δ%.1f°, Elong=%.2f | Funnel: Δ%.1f°, Elong=%.2f', ...
            frame, time_ms, mantle_rotation_angle(frame), mantle_elong, ...
            funnel_rotation_angle(frame), funnel_elong), 'FontSize', 9, 'Color', 'white');

        % Capture frame
        drawnow;
        frame_data = getframe(fig);
        frame_img = frame_data.cdata;

        % Initialize video writer on first frame
        if idx == 1
            [frame_h, frame_w, ~] = size(frame_img);
            video_writer = VideoWriter(video_filename, 'MPEG-4');
            video_writer.FrameRate = 10;
            open(video_writer);
            fprintf('Video writer initialized: %d x %d\n', frame_w, frame_h);
            expected_h = frame_h;
            expected_w = frame_w;
        end

        % Ensure frame matches expected dimensions
        [current_h, current_w, ~] = size(frame_img);
        if current_h ~= expected_h || current_w ~= expected_w
            frame_img = imresize(frame_img, [expected_h, expected_w]);
        end

        % Write to video
        writeVideo(video_writer, frame_img);

        % Save as PNG
        output_filename = sprintf('pca_frame_%04d.png', frame);
        output_path = fullfile(frames_dir, output_filename);
        imwrite(frame_img, output_path);

        close(fig);

        if mod(idx, 10) == 0 || idx == length(frames_to_render)
            fprintf('  Rendered %d/%d frames\n', idx, length(frames_to_render));
        end
    end

    % Close video writer
    if ~isempty(video_writer)
        close(video_writer);
    end

    fprintf('Frames saved to: %s\n', frames_dir);
    fprintf('Video created: %s\n', video_filename);
end

function [avg_dist_x, avg_dist_y, rotation_angle] = compute_rotated_lab_distances(tracks, pca_coeff1, pca_center, num_frames)
    % Compute average pairwise distances in rotated lab frame
    % Lab starts aligned with image axes at frame 1
    % Then rotates by the relative change in PCA angle

    num_tracks = length(tracks);

    % Get initial PCA angle (frame 1)
    initial_angle = atan2(pca_coeff1(1, 2), pca_coeff1(1, 1));

    % Preallocate output arrays
    avg_dist_x = zeros(num_frames, 1);
    avg_dist_y = zeros(num_frames, 1);
    rotation_angle = zeros(num_frames, 1);

    for frame = 1:num_frames
        % Get positions at this frame
        positions = [];
        for i = 1:num_tracks
            track = tracks{i};
            if frame <= size(track, 1)
                positions = [positions; track(frame, 2:3)];
            end
        end

        % Need at least 2 points for distance calculation
        if size(positions, 1) < 2
            continue;
        end

        % Get current PCA angle
        current_angle = atan2(pca_coeff1(frame, 2), pca_coeff1(frame, 1));

        % Relative rotation from initial frame
        relative_angle = current_angle - initial_angle;
        rotation_angle(frame) = rad2deg(relative_angle);

        % Get center position
        center = pca_center(frame, :);

        % Translate positions to center origin
        centered_positions = positions - center;

        % Rotate by NEGATIVE relative angle to align with lab frame
        % Lab frame rotates BY the relative angle, so we rotate positions back
        theta = -relative_angle;
        R = [cos(theta), -sin(theta);
             sin(theta),  cos(theta)];

        % Positions in lab frame (rotated from image axes by relative_angle)
        lab_positions = (R * centered_positions')';

        % Compute all pairwise distances along each axis
        n = size(lab_positions, 1);
        
        if n >= 2
            % Compute pairwise distances
            dist_x_all = [];
            dist_y_all = [];
            
            for i = 1:n
                for j = i+1:n
                    % Distance along lab X axis (horizontal)
                    dist_x_all = [dist_x_all; abs(lab_positions(i, 1) - lab_positions(j, 1))];

                    % Distance along lab Y axis (vertical)
                    dist_y_all = [dist_y_all; abs(lab_positions(i, 2) - lab_positions(j, 2))];
                end
            end
            
            % Compute average distances
            avg_dist_x(frame) = mean(dist_x_all);
            avg_dist_y(frame) = mean(dist_y_all);
        end
        
        if mod(frame, 50) == 0 || frame == num_frames
            fprintf('  Processed frame %d/%d\n', frame, num_frames);
        end
    end
    
    fprintf('  Rotation analysis complete!\n');
end
