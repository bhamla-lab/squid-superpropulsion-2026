% C20: Load frames, detect green boundary, and visualize
% 각 프레임 불러와서 초록색 검출 및 가장자리 계산, 가시화

clc; clear all; close all;

%% Parameters
VIDEO_NAME = 'S_F05mm_C01';  % Must match C01
MAT_FILE = fullfile('../03_output/00_start_detection', sprintf('%s.mat', VIDEO_NAME));
SNAPSHOT_FOLDER = fullfile('../03_output/04_snapshots', VIDEO_NAME);
OUTPUT_FOLDER = fullfile('../03_output/06_boundary_visualization', VIDEO_NAME);

% Green detection parameters (HSV)
GREEN_H_RANGE = [0.2, 0.5];   % Hue range
GREEN_S_MIN = 0.2;             % Saturation minimum
GREEN_V_MIN = 0.2;             % Value minimum

%% Load start info
if ~exist(MAT_FILE, 'file')
    error('Mat file not found: %s. Run C01 first.', MAT_FILE);
end
load(MAT_FILE, 'startInfo');

%% Setup
if ~exist(SNAPSHOT_FOLDER, 'dir')
    error('Snapshot folder not found: %s. Run C10 first.', SNAPSHOT_FOLDER);
end

if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

imageFiles = dir(fullfile(SNAPSHOT_FOLDER, '*.png'));

if isempty(imageFiles)
    error('No images found in: %s', SNAPSHOT_FOLDER);
end

fprintf('Processing %d images from: %s\n', length(imageFiles), SNAPSHOT_FOLDER);
fprintf('Saving results to: %s\n\n', OUTPUT_FOLDER);

%% Process each frame
for i = 1:length(imageFiles)
    % Load image
    filename = imageFiles(i).name;
    filepath = fullfile(SNAPSHOT_FOLDER, filename);
    img = imread(filepath);

    % Extract time from filename (e.g., S_F05mm_C01_t3s.png -> 3)
    tokens = regexp(filename, '_t(\d+)s\.png', 'tokens');
    if ~isempty(tokens)
        timeFromGreen = str2double(tokens{1}{1});
        actualTime = startInfo.greenTime + timeFromGreen;
    else
        actualTime = NaN;
    end

    % Detect green pixels
    [greenMask, greenBoundary] = detectGreenBoundary(img, GREEN_H_RANGE, GREEN_S_MIN, GREEN_V_MIN);

    % Visualize
    fig = figure('Position', [100, 100, 1400, 500], 'Visible', 'off');

    % Original image
    subplot(1, 3, 1);
    imshow(img);
    title('Original Frame', 'FontSize', 12);

    % Green mask
    subplot(1, 3, 2);
    imshow(greenMask);
    title('Green Detection Mask', 'FontSize', 12);

    % Overlay boundary on original
    subplot(1, 3, 3);
    imshow(img); hold on;
    if ~isempty(greenBoundary)
        plot(greenBoundary(:,2), greenBoundary(:,1), 'r-', 'LineWidth', 2);
    end
    title('Green Boundary Overlay', 'FontSize', 12);

    % Title with time info
    if ~isnan(actualTime)
        sgtitle(sprintf('%s | t=%.2fs (%.0fs after green)', filename, actualTime, timeFromGreen), ...
                'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
    else
        sgtitle(filename, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
    end

    % Save figure
    [~, baseName, ~] = fileparts(filename);
    outputFile = fullfile(OUTPUT_FOLDER, [baseName '_boundary.png']);
    saveas(fig, outputFile);
    close(fig);

    fprintf('Processed: %s (t=%.2fs) -> saved to 06_boundary_visualization\n', filename, actualTime);
end

fprintf('\nDone! Processed %d frames\n', length(imageFiles));

%% Helper Function
function [greenMask, boundary] = detectGreenBoundary(img, hRange, sMin, vMin)
    % Convert to HSV
    hsv = rgb2hsv(img);

    % Create green mask
    greenMask = (hsv(:,:,1) >= hRange(1) & hsv(:,:,1) <= hRange(2)) & ...
                (hsv(:,:,2) >= sMin) & ...
                (hsv(:,:,3) >= vMin);

    % Find boundary
    if any(greenMask(:))
        % Clean up mask
        greenMask = bwareaopen(greenMask, 50);  % Remove small regions
        greenMask = imfill(greenMask, 'holes'); % Fill holes

        % Get all boundaries
        boundaries = bwboundaries(greenMask, 'noholes');

        if ~isempty(boundaries)
            % Find region with largest width (주요영역)
            maxWidth = 0;
            maxIdx = 1;

            for i = 1:length(boundaries)
                b = boundaries{i};
                width = max(b(:,2)) - min(b(:,2));  % x range
                if width > maxWidth
                    maxWidth = width;
                    maxIdx = i;
                end
            end

            boundary = boundaries{maxIdx};  % Take widest region
        else
            boundary = [];
        end
    else
        boundary = [];
    end
end
