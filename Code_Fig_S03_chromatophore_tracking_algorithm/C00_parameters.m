function p = C00_parameters(dataset_index)
% C00_parameters.m
% Centralized parameters for squid deformation image processing
% Automatically detects and selects datasets from 01_raw_data/ folder
% Author: Claude
% Date: 2025-10-13
%
% Usage:
%   p = C00_parameters(1);     % Select first dataset
%   p = C00_parameters(2);     % Select second dataset
%   p = C00_parameters();      % Use default (first) dataset
%
% The function automatically:
%   - Scans 01_raw_data/ for available datasets
%   - Detects file patterns and extensions
%   - Finds first image file
%   - Sets up all output directories

% ========== SCAN AVAILABLE DATASETS ==========
raw_data_dir = '01_raw_data';

% Get all folders in raw_data directory
if ~exist(raw_data_dir, 'dir')
    error('Raw data directory not found: %s\nPlease create this directory and add your datasets.', raw_data_dir);
end

all_items = dir(raw_data_dir);

% Separate folders and video files
dataset_folders = all_items([all_items.isdir] & ~ismember({all_items.name}, {'.', '..', '.DS_Store'}));
video_files = all_items(~[all_items.isdir] & endsWith({all_items.name}, {'.mp4', '.avi', '.mov'}));

% Combine: folders first, then video files
all_datasets = [dataset_folders; video_files];

if isempty(all_datasets)
    error('No datasets or videos found in %s\nPlease add dataset folders or video files to this directory.', raw_data_dir);
end

% Sort alphabetically
[~, sort_idx] = sort({all_datasets.name});
all_datasets = all_datasets(sort_idx);

% ========== SELECT DATASET ==========
if nargin < 1 || isempty(dataset_index)
    dataset_index = 1;  % Default to first dataset
end

if dataset_index < 1 || dataset_index > length(all_datasets)
    error('Invalid dataset index: %d\nAvailable datasets: 1-%d', dataset_index, length(all_datasets));
end

selected_item = all_datasets(dataset_index);
selected_dataset = selected_item.name;
is_video = ~selected_item.isdir;

% ========== DISPLAY AVAILABLE DATASETS ==========
fprintf('\n========== DATASET SELECTION ==========\n');
fprintf('Available datasets in %s/:\n', raw_data_dir);

for i = 1:length(all_datasets)
    item = all_datasets(i);
    item_name = item.name;

    if item.isdir
        % It's a folder - count image files
        dataset_path = fullfile(raw_data_dir, item_name);
        files = dir(fullfile(dataset_path, '*.*'));
        image_files = files(~[files.isdir] & ~startsWith({files.name}, '.'));

        if ~isempty(image_files)
            extensions = cellfun(@(x) extractAfter(x, '.'), {image_files.name}, 'UniformOutput', false);
            [unique_ext, ~, idx] = unique(extensions);
            counts = accumarray(idx(:), 1);
            [max_count, max_idx] = max(counts);
            main_ext = unique_ext{max_idx};

            if i == dataset_index
                fprintf('  [%d] %s (%d .%s files) <- SELECTED\n', i, item_name, max_count, main_ext);
            else
                fprintf('  [%d] %s (%d .%s files)\n', i, item_name, max_count, main_ext);
            end
        else
            if i == dataset_index
                fprintf('  [%d] %s (no files found) <- SELECTED\n', i, item_name);
            else
                fprintf('  [%d] %s (no files found)\n', i, item_name);
            end
        end
    else
        % It's a video file
        if i == dataset_index
            fprintf('  [%d] %s (VIDEO) <- SELECTED\n', i, item_name);
        else
            fprintf('  [%d] %s (VIDEO)\n', i, item_name);
        end
    end
end

fprintf('\nSelected dataset: %s\n', selected_dataset);

% ========== AUTO-DETECT FILE PATTERN OR VIDEO INFO ==========
if is_video
    % Handle video file
    video_path = fullfile(raw_data_dir, selected_dataset);
    video_obj = VideoReader(video_path);

    fprintf('Video file detected\n');
    fprintf('Duration: %.2f seconds\n', video_obj.Duration);
    fprintf('Frame rate: %.2f fps\n', video_obj.FrameRate);
    fprintf('Total frames: %d\n', round(video_obj.Duration * video_obj.FrameRate));
    fprintf('Resolution: %dx%d\n', video_obj.Width, video_obj.Height);

    % Set video-specific fields
    image_pattern = '';
    first_image = '';
    detected_extension = 'video';
    num_images = round(video_obj.Duration * video_obj.FrameRate);

else
    % Handle image folder
    dataset_path = fullfile(raw_data_dir, selected_dataset);

    % Find all image files
    all_files = dir(dataset_path);
    all_files = all_files(~[all_files.isdir] & ~startsWith({all_files.name}, '.'));

    if isempty(all_files)
        error('No files found in dataset folder: %s', dataset_path);
    end

    % Get file extensions
    file_extensions = cellfun(@(x) extractAfter(x, '.'), {all_files.name}, 'UniformOutput', false);

    % Find most common extension
    [unique_ext, ~, idx] = unique(file_extensions);
    counts = accumarray(idx(:), 1);
    [~, max_idx] = max(counts);
    detected_extension = unique_ext{max_idx};

    % Filter files by most common extension
    image_files = all_files(strcmp(file_extensions, detected_extension));

    if isempty(image_files)
        error('No image files with extension .%s found in: %s', detected_extension, dataset_path);
    end

    % Sort files by name
    [~, sort_idx] = sort({image_files.name});
    image_files = image_files(sort_idx);

    % Get first image
    first_image = image_files(1).name;

    % Extract base pattern (remove sequential numbers and extension)
    % Try to detect pattern by finding common prefix
    all_names = {image_files.name};
    if length(all_names) > 1
        % Find common prefix
        first_name = all_names{1};
        last_name = all_names{end};

        common_len = 0;
        for k = 1:min(length(first_name), length(last_name))
            if first_name(k) == last_name(k)
                common_len = k;
            else
                break;
            end
        end

        % Pattern is common prefix + wildcard + extension
        if common_len > 0
            base_pattern = first_name(1:common_len);
            image_pattern = [base_pattern '*.' detected_extension];
        else
            % Fallback: use extension only
            image_pattern = ['*.' detected_extension];
        end
    else
        % Only one file: use extension
        image_pattern = ['*.' detected_extension];
    end

    fprintf('File pattern: %s\n', image_pattern);
    fprintf('First image: %s\n', first_image);
    fprintf('Total images: %d\n', length(image_files));

    num_images = length(image_files);
end

% ========== SET DATA PATHS ==========
p.dataset_name = selected_dataset;
p.dataset_index = dataset_index;
p.is_video = is_video;

if is_video
    p.data_dir = raw_data_dir;
    p.video_path = fullfile(raw_data_dir, selected_dataset);
    [~, dataset_basename, ~] = fileparts(selected_dataset);
    p.output_dir_02 = fullfile('02_processed_data', dataset_basename);
    p.output_dir_03 = fullfile('03_cleaned_data', dataset_basename);
    p.output_dir_04 = fullfile('04_tracking_data', dataset_basename);
    p.output_dir_05 = fullfile('05_combined_analysis', dataset_basename);
    p.output_dir_16 = fullfile('16_pca_rotation_analysis', dataset_basename);
else
    p.data_dir = fullfile(raw_data_dir, selected_dataset);
    p.output_dir_02 = fullfile('02_processed_data', selected_dataset);
    p.output_dir_03 = fullfile('03_cleaned_data', selected_dataset);
    p.output_dir_04 = fullfile('04_tracking_data', selected_dataset);
    p.output_dir_05 = fullfile('05_combined_analysis', selected_dataset);
    p.output_dir_16 = fullfile('16_pca_rotation_analysis', selected_dataset);
end

p.image_pattern = image_pattern;
p.first_image_name = first_image;
p.image_extension = detected_extension;
p.num_images = num_images;

% ========== DATASET-SPECIFIC PARAMETERS ==========
% Configure different parameters for each dataset
% Add new datasets here with their custom parameters

fprintf('\n========== LOADING PARAMETERS ==========\n');

% Define default parameters (used if dataset not found below)
default_params = struct();

% C01, C02 parameters
default_params.gaussian_sigma = 30;
default_params.threshold_percentile = 2;
default_params.min_dot_area = 1;
default_params.max_dot_area = 100;
default_params.save_format = 'png';
default_params.preview_mode = false;
default_params.show_progress_every = 50;

% C03 parameters
default_params.num_frames_to_check = 100;
default_params.stationary_threshold = 0.8;
default_params.flicker_threshold = 0.2;
default_params.flicker_max_size = 2;

% C04 parameters
default_params.max_displacement = 5;
default_params.min_track_length_c04 = 5;
default_params.num_frames_to_process = 0;

% C05 parameters
default_params.max_displacement_per_frame = 5;
default_params.min_track_length_c05 = 300;
default_params.min_total_displacement = 5;

% C06 parameters
default_params.frame_rate = 1000;
default_params.require_full_duration_tracks = true;
default_params.num_frames_to_analyze = 0;  % 0 = all frames, n = process first n frames in C16

% C08 parameters (time normalization)
default_params.time_normalization = 0;  % Time offset in ms to normalize all datasets

% Dataset-specific parameter overrides
switch selected_dataset
    case 'Dataset_00_10_0000_0000'  % Index: 1
        fprintf('Using parameters for: %s (Index: 1)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 0;

    case 'Dataset_09_02_6613_7688'  % Index: 2
        fprintf('Using parameters for: %s (Index: 2)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 103;

    case 'Dataset_09_03_0358_0749'  % Index: 3
        fprintf('Using parameters for: %s (Index: 3)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 0;

    case 'Dataset_09_03_11203_12147'  % Index: 4
        fprintf('Using parameters for: %s (Index: 4)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 295;

    case 'Dataset_09_03_5340_6089'  % Index: 5
        fprintf('Using parameters for: %s (Index: 5)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 0;

    case 'Dataset_09_04_0838_1416'  % Index: 6
        fprintf('Using parameters for: %s (Index: 6)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 0;

    case 'Dataset_10_01_0521_1172'  % Index: 7
        fprintf('Using parameters for: %s (Index: 7)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 46;

    case 'Dataset_10_02_0228_0814'  % Index: 8
        fprintf('Using parameters for: %s (Index: 8)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 91;

    case 'Dataset_10_03_0391_0830'  % Index: 9 (Original reference dataset)
        fprintf('Using parameters for: %s (Index: 9, Reference)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 0;

    case 'Dataset_12_01_0894_1448'  % Index: 29
        fprintf('Using parameters for: %s (Index: 29)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 120;

    case 'Dataset_12_01_7340_8023'  % Index: 30
        fprintf('Using parameters for: %s (Index: 30)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 0;

    case 'Dataset_12_02_0364_0610'  % Index: 31
        fprintf('Using parameters for: %s (Index: 31)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 200;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 31;

    case 'Dataset_12_03_7232_7356'  % Index: 32
        fprintf('Using parameters for: %s (Index: 32)\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 1;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 70;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 0;

    case 'Dataset_12_07_1116_1995_S'  % Index: 33 (Severed case)
        fprintf('Using parameters for: %s (Index: 33) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_12_07_3199_3948_S'  % Index: 34 (Severed case)
        fprintf('Using parameters for: %s (Index: 34) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_12_07_4306_5935_S'  % Index: 35 (Severed case)
        fprintf('Using parameters for: %s (Index: 35) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_12_08_0000_1140_S'  % Index: 36 (Severed case)
        fprintf('Using parameters for: %s (Index: 36) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_12_08_1108_1498_S'  % Index: 37 (Severed case)
        fprintf('Using parameters for: %s (Index: 37) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 0;

    case 'Dataset_12_08_2019_2768_S'  % Index: 38 (Severed case)
        fprintf('Using parameters for: %s (Index: 38) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 123;


    case 'Dataset_12_08_3190_3842_S'  % Index: 39 (Severed case)
        fprintf('Using parameters for: %s (Index: 39) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 194;

    case 'Dataset_12_08_4590_5437_S'  % Index: 40 (Severed case)
        fprintf('Using parameters for: %s (Index: 40) - SEVERED\n', selected_dataset);
        params = default_params;
        params.time_normalization = 266;


    case 'Dataset_10_07_0163_1270_S'  % Index: 10 (Severed case)
        fprintf('Using parameters for: %s (Index: 10) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_10_07_2011_2924_S'  % Index: 11 (Severed case)
        fprintf('Using parameters for: %s (Index: 11) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 2;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_10_07_3345_4290_S'  % Index: 12 (Severed case)
        fprintf('Using parameters for: %s (Index: 12) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 308;
    case 'Dataset_10_07_4811_5788_S'  % Index: 13 (Severed case)
        fprintf('Using parameters for: %s (Index: 13) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 300;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_10_07_5893_7156_S'  % Index: 14 (Severed case)
        fprintf('Using parameters for: %s (Index: 14) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_10_07_7612_8947_S'  % Index: 15 (Severed case)
        fprintf('Using parameters for: %s (Index: 15) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_10_08_0000_1238_S'  % Index: 16 (Severed case)
        fprintf('Using parameters for: %s (Index: 16) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 200;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_10_08_1336_2378_S'  % Index: 17 (Severed case)
        fprintf('Using parameters for: %s (Index: 17) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 200;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_10_08_2801_3843_S'  % Index: 18 (Severed case)
        fprintf('Using parameters for: %s (Index: 18) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 200;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_10_08_3515_4396_S'  % Index: 19 (Severed case)
        fprintf('Using parameters for: %s (Index: 19) - SEVERED\n', selected_dataset);
       params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 200;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_10_08_4883_6447_S'  % Index: 20 (Severed case)
        fprintf('Using parameters for: %s (Index: 20) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 200;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_10_08_6739_7911_S'  % Index: 21 (Severed case)
        fprintf('Using parameters for: %s (Index: 21) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_11_06_0358_1325_S'  % Index: 22 (Severed case)
        fprintf('Using parameters for: %s (Index: 22) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 5;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 5;
        params.min_track_length_c05 = 200;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_11_06_1487_2302_S'  % Index: 23 (Severed case)
        fprintf('Using parameters for: %s (Index: 23) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 55;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 800;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_11_06_2855_3636_S'  % Index: 24 (Severed case)
        fprintf('Using parameters for: %s (Index: 24) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 55;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 350;  % 0 = all frames, n = process first n frames in C16
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_11_06_4092_4613_S'  % Index: 25 (Severed case)
        fprintf('Using parameters for: %s (Index: 25) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 55;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 0;  % 0 = all frames, n = process first n frames in C16
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_11_06_5231_6078_S'  % Index: 26 (Severed case)
        fprintf('Using parameters for: %s (Index: 26) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 15;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 15;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 0;  % 0 = all frames, n = process first n frames in C16
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_11_06_7706_8944_S'  % Index: 27 (Severed case)
        fprintf('Using parameters for: %s (Index: 27) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 500;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 15;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 15;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 0;  % 0 = all frames, n = process first n frames in C16
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_11_06_a735_b582_S'  % Index: 28 (Severed case)
        fprintf('Using parameters for: %s (Index: 28) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_08_04_0554_1498'  % Index: 2
        fprintf('Using parameters for: %s (Index: 2)\n', selected_dataset);
        params = default_params;
        % C01, C02 parameters
        params.gaussian_sigma = 10;
        % C03 parameters
        % params.num_frames_to_check = 100;
        params.stationary_threshold = 0.9;
        % params.flicker_threshold = 0.2;
        % params.flicker_max_size = 2;
        % % C04 parameters
        % params.max_displacement = 5;
        % params.min_track_length_c04 = 5;
        % params.num_frames_to_process = 0;
        % % C05 parameters
        % params.max_displacement_per_frame = 5;
        % params.min_track_length_c05 = 70;
        % params.min_total_displacement = 5;
        % % C06 parameters
        % params.frame_rate = 1000;
        % params.require_full_duration_tracks = true;
        % C08 parameters
        params.time_normalization = 0;
    case 'Dataset_08_04_1693_2866'  % Index: 3
        fprintf('Using parameters for: %s (Index: 3)\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
    case 'Dataset_08_05_2867_4007'  % Index: 4
        fprintf('Using parameters for: %s (Index: 4)\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
    case 'Dataset_08_06_0000_1173'  % Index: 5
        fprintf('Using parameters for: %s (Index: 5)\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
    case 'Dataset_08_06_1205_2313'  % Index: 6
        fprintf('Using parameters for: %s (Index: 6)\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
    case 'Dataset_08_06_3387_5081'  % Index: 7
        fprintf('Using parameters for: %s (Index: 7)\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
    case 'Dataset_08_08_0293_1563_S'  % Index: 8 (Severed case)
        fprintf('Using parameters for: %s (Index: 8) - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
    case 'Dataset_08_08_1530_2638_S'  % Index: 9 (Severed case)
        fprintf('Using parameters for: %s (Index: 9) - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
    case 'Dataset_08_08_2572_3451_S'  % Index: 10 (Severed case)
        fprintf('Using parameters for: %s (Index: 10) - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
    case 'Dataset_08_08_3809_5014_S'  % Index: 11 (Severed case)
        fprintf('Using parameters for: %s (Index: 11) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_08_09_0000_0847_S'  % Index: 12 (Severed case - Dislocation Funnel)
        fprintf('Using parameters for: %s (Index: 12) - SEVERED (Dislocation Funnel)\n', selected_dataset);
        params = default_params;

    case 'Dataset_08_09_10608_12171_S'  % Index: 13 (Severed case)
        fprintf('Using parameters for: %s (Index: 13) - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 300;  % 0 = all frames, n = process first n frames in C16

    case 'Dataset_08_09_1108_1922_S'  % Index: 14 (Severed case)
        fprintf('Using parameters for: %s (Index: 14) - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
    case 'Dataset_08_09_2524_3292_S'  % Index: 15 (Severed case)
        fprintf('Using parameters for: %s (Index: 15) - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 0;
        params.num_frames_to_analyze = 600;  % 0 = all frames, n = process first n frames in C16

    case 'Dataset_08_09_3705_4421_S'  % Index: 16 (Severed case)
        fprintf('Using parameters for: %s (Index: 16) - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 600;  % 0 = all frames, n = process first n frames in C16
    case 'Dataset_08_09_4715_5431_S'  % Index: 17 (Severed case)
        fprintf('Using parameters for: %s (Index: 17) - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 600;  % 0 = all frames, n = process first n frames in C16
    case 'Dataset_08_09_5463_6766_S'  % Index: 18 (Severed case)
        fprintf('Using parameters for: %s (Index: 18) - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 600;  % 0 = all frames, n = process first n frames in C16
    case 'Dataset_08_09_6798_8133_S'  % Index: 19 (Severed case)
        fprintf('Using parameters for: %s (Index: 19) - SEVERED\n', selected_dataset);
        params = default_params;
    case 'Dataset_08_09_8003_9403_S'  % Index: 20 (Severed case)
        fprintf('Using parameters for: %s (Index: 20) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_08_09_9468_10543_S'  % Index: 21 (Severed case)
        fprintf('Using parameters for: %s (Index: 21) - SEVERED\n', selected_dataset);
        params = default_params;

    case 'Dataset_08_10_1051_1503_S'  % (Severed case)
        fprintf('Using parameters for: %s - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 600;  % 0 = all frames, n = process first n frames in C16

    case 'Dataset_08_10_1702_2277_S'  % (Severed case)
        fprintf('Using parameters for: %s - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 600;  % 0 = all frames, n = process first n frames in C16

    case 'Dataset_08_10_2636_3384_S'  % (Severed case)
        fprintf('Using parameters for: %s - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 600;  % 0 = all frames, n = process first n frames in C16

    case 'Dataset_08_10_3624_4243_S'  % (Severed case)
        fprintf('Using parameters for: %s - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 600;  % 0 = all frames, n = process first n frames in C16
    case 'Dataset_08_10_4675_5684_S'  % (Severed case)
        fprintf('Using parameters for: %s - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 600;  % 0 = all frames, n = process first n frames in C16
    case 'squid08'  % Index: 41 (Video frames folder)
        fprintf('Using parameters for: %s (Index: 41) - VIDEO FRAMES\n', selected_dataset);
        params = default_params;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.8;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 50;
        params.min_track_length_c05 = 120;
        params.min_total_displacement = 5;
        % C06 parameters
        params.require_full_duration_tracks = false;
        params.time_normalization = 0;
        params.frame_rate = 30;


    case 'squid08.mp4'  % Index: 42 (Video file)
        fprintf('Using parameters for: %s (Index: 41) - VIDEO\n', selected_dataset);
        params = default_params;
        params.min_track_length_c05 = 10;
        params.require_full_duration_tracks = false;
        params.time_normalization = 0;


    case 'squid09'  % Index: 43 (Video frames folder)
        fprintf('Using parameters for: %s (Index: 43) - VIDEO FRAMES\n', selected_dataset);
        params = default_params;
        params.frame_rate = 30;

    case 'squid09.mp4'  % Index: 44 (Video file)
        fprintf('Using parameters for: %s (Index: 44) - VIDEO\n', selected_dataset);
        params = default_params;
        params.frame_rate = 30;

    case 'squid10'  % Index: 45 (Video frames folder)
        fprintf('Using parameters for: %s (Index: 45) - VIDEO FRAMES\n', selected_dataset);
        params = default_params;
        params.frame_rate = 30;

    case 'squid10.mp4'  % Index: 46 (Video file)
        fprintf('Using parameters for: %s (Index: 46) - VIDEO\n', selected_dataset);
        params = default_params;
        params.frame_rate = 30;

    case 'squid11'  % Index: 47 (Video frames folder)
        fprintf('Using parameters for: %s (Index: 47) - VIDEO FRAMES\n', selected_dataset);
        params = default_params;
        params.frame_rate = 30;

    case 'squid11.mp4'  % Index: 48 (Video file)
        fprintf('Using parameters for: %s (Index: 48) - VIDEO\n', selected_dataset);
        params = default_params;
        params.frame_rate = 30;

    % ========== SQUID 3 (Dataset_03_*) ==========
    case 'Dataset_03_01_0252_1120'  % Index: 2
        fprintf('Using parameters for: %s (Index: 2) - Squid 3\n', selected_dataset);
        params = default_params;
        % C01, C02 parameters
        % params.gaussian_sigma = 30;
        params.threshold_percentile = 10;
        % params.min_dot_area = 1;
        % params.max_dot_area = 100;
        % params.save_format = 'png';
        % params.preview_mode = false;
        % params.show_progress_every = 50;
        % % C03 parameters
        % params.num_frames_to_check = 100;
        % params.stationary_threshold = 0.8;
        % params.flicker_threshold = 0.2;
        % params.flicker_max_size = 2;
        % % C04 parameters
        % params.max_displacement = 2;
        % params.min_track_length_c04 = 5;
        % params.num_frames_to_process = 0;
        % % C05 parameters
        % params.max_displacement_per_frame = 5;
        % params.min_track_length_c05 = 500;
        % params.min_total_displacement = 5;
        % % C06 parameters
        % params.frame_rate = 1000;
        % params.require_full_duration_tracks = false;
        % % C08 parameters
        % params.time_normalization = 0;

    case 'Dataset_03_01_2099_3247'  % Index: 3
        fprintf('Using parameters for: %s (Index: 3) - Squid 3\n', selected_dataset);
        params = default_params;
        params.threshold_percentile = 10;

    case 'Dataset_03_02_0364_1512'  % Index: 4
        fprintf('Using parameters for: %s (Index: 4) - Squid 3\n', selected_dataset);
        params = default_params;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_03_0000_1036'  % Index: 5
        fprintf('Using parameters for: %s (Index: 5) - Squid 3\n', selected_dataset);
        params = default_params;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_03_0924_1708'  % Index: 6
        fprintf('Using parameters for: %s (Index: 6) - Squid 3\n', selected_dataset);
        params = default_params;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_03_2855_3499'  % Index: 7
        fprintf('Using parameters for: %s (Index: 7) - Squid 3\n', selected_dataset);
        params = default_params;
        params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_03_3835_5067'  % Index: 8
        fprintf('Using parameters for: %s (Index: 8) - Squid 3\n', selected_dataset);
        params = default_params;
                params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_03_6186_7250'  % Index: 9
        fprintf('Using parameters for: %s (Index: 9) - Squid 3\n', selected_dataset);
        params = default_params;
                params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_03_9881_10861'  % Index: 10
        fprintf('Using parameters for: %s (Index: 10) - Squid 3\n', selected_dataset);
        params = default_params;
                params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_04_0252_0924'  % Index: 11
        fprintf('Using parameters for: %s (Index: 11) - Squid 3\n', selected_dataset);
        params = default_params;
                params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_04_0812_1524'  % Index: 12
        fprintf('Using parameters for: %s (Index: 12) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;
    case 'Dataset_03_04_2408_3136'  % Index: 13
        fprintf('Using parameters for: %s (Index: 13) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;
    case 'Dataset_03_04_3164_3975'  % Index: 14
        fprintf('Using parameters for: %s (Index: 14) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;
    case 'Dataset_03_04_3975_4899'  % Index: 15
        fprintf('Using parameters for: %s (Index: 15) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_04_4927_6186'  % Index: 16
        fprintf('Using parameters for: %s (Index: 16) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_04_6717_7921'  % Index: 17
        fprintf('Using parameters for: %s (Index: 17) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 10;
        params.max_displacement_per_frame = 2;
        
    case 'Dataset_03_05_0000_0868'  % Index: 18
        fprintf('Using parameters for: %s (Index: 18) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;
        
    case 'Dataset_03_05_0812_1708'  % Index: 19
        fprintf('Using parameters for: %s (Index: 19) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;
        
    case 'Dataset_03_05_1708_2576'  % Index: 20
        fprintf('Using parameters for: %s (Index: 20) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;
        
    case 'Dataset_03_05_2576_3471'  % Index: 21
        fprintf('Using parameters for: %s (Index: 21) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;
        
    case 'Dataset_03_05_3946_5009'  % Index: 22
        fprintf('Using parameters for: %s (Index: 22) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_05_8648_9544'  % Index: 23
        fprintf('Using parameters for: %s (Index: 23) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_05_9544_10917'  % Index: 24
        fprintf('Using parameters for: %s (Index: 24) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_06_0280_1036'  % Index: 25
        fprintf('Using parameters for: %s (Index: 25) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_06_1036_1624'  % Index: 26
        fprintf('Using parameters for: %s (Index: 26) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_06_2072_3052'  % Index: 27
        fprintf('Using parameters for: %s (Index: 27) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_06_8426_9256'  % Index: 28
        fprintf('Using parameters for: %s (Index: 28) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_06_9628_10692'  % Index: 29
        fprintf('Using parameters for: %s (Index: 29) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_07_0663_1482'  % Index: 30
        fprintf('Using parameters for: %s (Index: 30) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_07_1521_2339'  % Index: 31
        fprintf('Using parameters for: %s (Index: 31) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_07_2417_3274'  % Index: 32
        fprintf('Using parameters for: %s (Index: 32) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_07_5569_6545'  % Index: 33
        fprintf('Using parameters for: %s (Index: 33) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_07_3389_4754'  % Index: 34
        fprintf('Using parameters for: %s (Index: 34) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_07_9155_9817'  % Index: 35
        fprintf('Using parameters for: %s (Index: 35) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_08_1560_2261'  % Index: 36
        fprintf('Using parameters for: %s (Index: 36) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_08_2572_3391'  % Index: 37
        fprintf('Using parameters for: %s (Index: 37) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_03_08_4130_5066'  % Index: 38
        fprintf('Using parameters for: %s (Index: 38) - Squid 3\n', selected_dataset);
        params = default_params;
          params.frame_rate = 1000;
        params.threshold_percentile = 8;
        params.max_displacement_per_frame = 2;

    case 'Dataset_09_04_0838_1416_sideview'  % Index: 39
        fprintf('Using parameters for: %s (Index: 39) - Sideview\n', selected_dataset);
        params = default_params;
        params.frame_rate = 1000;
    case 'Dataset_11_05_1142_2040_S'  % (Severed case)
        fprintf('Using parameters for: %s - SEVERED\n', selected_dataset);
        params = default_params;
        params.gaussian_sigma = 10;
        params.stationary_threshold = 0.9;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 600;  % 0 = all frames, n = process first n frames in C16

    case 'Dataset_02_12_0728_1540'  % Index: 23 (Severed case)
        fprintf('Using parameters for: %s (Index: 23) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 55;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 800;
        % C08 parameters
        params.time_normalization = 308;

        case 'Dataset_02_12_1651_2490'  % Index: 23 (Severed case)
        fprintf('Using parameters for: %s (Index: 23) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 55;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 800;
        % C08 parameters
        params.time_normalization = 308;

        case 'Dataset_02_12_3021_3665'  % Index: 23 (Severed case)
        fprintf('Using parameters for: %s (Index: 23) - SEVERED\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 100;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 55;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 800;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_02_14_0000_0728'
        fprintf('Using parameters for: %s\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 500;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 55;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 470;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_02_14_0812_1568'
        fprintf('Using parameters for: %s\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 500;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 55;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 800;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_02_14_1484_2015'
        fprintf('Using parameters for: %s\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 500;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 55;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 800;
        % C08 parameters
        params.time_normalization = 308;

    case 'Dataset_02_14_6576_7471'
        fprintf('Using parameters for: %s\n', selected_dataset);
        params = struct();
        % C01, C02 parameters
        params.gaussian_sigma = 30;
        params.threshold_percentile = 5;
        params.min_dot_area = 1;
        params.max_dot_area = 100;
        params.save_format = 'png';
        params.preview_mode = false;
        params.show_progress_every = 50;
        % C03 parameters
        params.num_frames_to_check = 500;
        params.stationary_threshold = 0.99;
        params.flicker_threshold = 0.2;
        params.flicker_max_size = 2;
        % C04 parameters
        params.max_displacement = 55;
        params.min_track_length_c04 = 5;
        params.num_frames_to_process = 0;
        % C05 parameters
        params.max_displacement_per_frame = 55;
        params.min_track_length_c05 = 100;
        params.min_total_displacement = 5;
        % C06 parameters
        params.frame_rate = 1000;
        params.require_full_duration_tracks = 1;
        params.num_frames_to_analyze = 800;
        % C08 parameters
        params.time_normalization = 308;

    otherwise
        % Unknown dataset - use defaults
        fprintf('Using DEFAULT parameters for: %s (Not configured yet)\n', selected_dataset);
        fprintf('  -> Add custom parameters in C00_parameters.m if needed\n');
        params = default_params;
end

% Apply parameters to output structure
% ========================================================================
% C01, C02: DOT DETECTION AND PREPROCESSING
% ========================================================================
p.gaussian_sigma = params.gaussian_sigma;          % Sigma for Gaussian blur (larger = smoother background)
                                                    % Typical range: 10-50 pixels

p.threshold_percentile = params.threshold_percentile;    % Detect pixels darker than this percentile (0-100)
                                                          % Lower = detect darker dots only
                                                          % Higher = detect more dots (including lighter ones)

p.min_dot_area = params.min_dot_area;             % Minimum area of dots to keep (in pixels)
                                                   % Removes noise and tiny artifacts

p.max_dot_area = params.max_dot_area;             % Maximum area of dots to keep (in pixels)
                                                   % Removes large blobs that aren't individual dots

p.save_format = params.save_format;          % Output format: 'png' or 'bmp'
                                             % PNG is compressed, BMP is uncompressed

p.preview_mode = params.preview_mode;         % true: process only first 10 images
                                              % false: process all images

p.show_progress_every = params.show_progress_every;     % Show progress every N frames during batch processing

% ========================================================================
% C03: STATIONARY NOISE DETECTION
% ========================================================================
p.num_frames_to_check = params.num_frames_to_check;      % Number of frames to analyze for noise (0 = all frames)
p.stationary_threshold = params.stationary_threshold;     % Fraction of frames where dot must appear (0-1)
                                                           % 0.8 means dot appears in 80%+ of frames = stationary
p.flicker_threshold = params.flicker_threshold;           % Minimum frequency to detect flickering dots (0-1)
                                                          % 0.2 means dot appears in 20%+ of frames
p.flicker_max_size = params.flicker_max_size;            % Maximum size (pixels) for flickering dot removal

% ========================================================================
% C04: DOT TRACKING
% ========================================================================
p.max_displacement = params.max_displacement;         % Maximum pixels a dot can move between frames
                                                       % Increase if dots move fast, decrease to reduce false matches

p.min_track_length_c04 = params.min_track_length_c04;   % Minimum frames to keep a track (initial filter)
                                                          % Short tracks are likely noise

p.num_frames_to_process = params.num_frames_to_process;  % Number of frames to track (0 = all frames)
                                                           % Useful for testing on subset

% ========================================================================
% C05: TRACK FILTERING AND REGION SELECTION
% ========================================================================
p.max_displacement_per_frame = params.max_displacement_per_frame;  % Max pixels allowed between consecutive positions
                                                                    % Rejects tracks with large jumps (likely mismatches)

p.min_track_length_c05 = params.min_track_length_c05;   % Minimum track length for analysis
                                                          % Should be >> min_track_length_c04

p.min_total_displacement = params.min_total_displacement;  % Minimum total displacement (pixels)
                                                            % Filters out stationary dots

% ========================================================================
% C06: DEFORMATION ANALYSIS
% ========================================================================
p.frame_rate = params.frame_rate;             % Frame rate in Hz (1000 Hz = 1 kHz = 1 ms per frame)
                                              % Used to convert frame numbers to time

p.require_full_duration_tracks = params.require_full_duration_tracks;  % true: Only use tracks spanning full duration
                                                                        % false: Use all tracks >= min_track_length_c05

p.num_frames_to_analyze = params.num_frames_to_analyze;  % 0 = all frames, n = process first n frames in C16

% ========================================================================
% C08: TIME NORMALIZATION
% ========================================================================
p.time_normalization = params.time_normalization;  % Time offset in ms to normalize all datasets

% Display loaded parameters
fprintf('\n========== PARAMETERS SUMMARY ==========\n');
fprintf('C01, C02 (Detection):\n');
fprintf('  Gaussian sigma: %d\n', p.gaussian_sigma);
fprintf('  Threshold percentile: %d\n', p.threshold_percentile);
fprintf('  Dot area: %d-%d pixels\n', p.min_dot_area, p.max_dot_area);
fprintf('  Save format: %s\n', p.save_format);
fprintf('\nC04 (Tracking):\n');
fprintf('  Max displacement: %d pixels\n', p.max_displacement);
fprintf('  Min track length (initial): %d frames\n', p.min_track_length_c04);
if p.num_frames_to_process > 0
    fprintf('  Frames to process: %d\n', p.num_frames_to_process);
else
    fprintf('  Frames to process: ALL\n');
end
fprintf('\nC05 (Filtering):\n');
fprintf('  Max displacement/frame: %d pixels\n', p.max_displacement_per_frame);
fprintf('  Min track length (final): %d frames\n', p.min_track_length_c05);
fprintf('  Min total displacement: %d pixels\n', p.min_total_displacement);
fprintf('\nC06 (Deformation):\n');
fprintf('  Frame rate: %d Hz\n', p.frame_rate);
if p.require_full_duration_tracks
    fprintf('  Track filter: Full-duration only\n');
else
    fprintf('  Track filter: All tracks >= %d frames\n', p.min_track_length_c05);
end
fprintf('========================================\n\n');

end
