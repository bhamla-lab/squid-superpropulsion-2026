% C10: Extract cropped snapshots at specified times after green detection
% 초록색 검출 시점부터 지정 시간 후 스냅샷 크롭하여 저장

clc; clear all; close all;

%% Parameters
VIDEO_NAME = 'S_F05mm_C01';  % Must match C01
MAT_FILE = fullfile('../03_output/00_start_detection', sprintf('%s.mat', VIDEO_NAME));
SNAPSHOT_TIMES = [0:2:20];  % Seconds after green detection
CROP_WIDTH_PERCENT = 38;   % Crop width as % of original image width
CROP_HEIGHT_PERCENT = 40;  % Crop height as % of original image height
OUTPUT_BASE = '../03_output/04_snapshots';

%% Load start info
if ~exist(MAT_FILE, 'file')
    error('Mat file not found: %s. Run C01 first.', MAT_FILE);
end
load(MAT_FILE, 'startInfo');

%% Setup output folder
outputFolder = fullfile(OUTPUT_BASE, VIDEO_NAME);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Open video
video = VideoReader(startInfo.videoPath);
fps = startInfo.fps;

fprintf('Video: %s\n', VIDEO_NAME);
fprintf('Green detected at: frame %d (%.2fs)\n', startInfo.greenFrame, startInfo.greenTime);
fprintf('Extracting snapshots at: %s seconds\n\n', mat2str(SNAPSHOT_TIMES));

%% Get crop parameters
% Read first frame to get dimensions
video.CurrentTime = 0;
tempFrame = readFrame(video);
[imgH, imgW, ~] = size(tempFrame);

% Calculate crop size
cropW = round(imgW * CROP_WIDTH_PERCENT / 100);
cropH = round(imgH * CROP_HEIGHT_PERCENT / 100);

% Center position from C01
cx = startInfo.greenPixelPos(1);
cy = startInfo.greenPixelPos(2);

fprintf('Image size: %d x %d\n', imgW, imgH);
fprintf('Crop size: %d x %d (%.0f%% x %.0f%%)\n', cropW, cropH, CROP_WIDTH_PERCENT, CROP_HEIGHT_PERCENT);
fprintf('Center: [%d, %d]\n\n', cx, cy);

%% Extract snapshots
for i = 1:length(SNAPSHOT_TIMES)
    dt = SNAPSHOT_TIMES(i);
    targetTime = startInfo.greenTime + dt;
    targetFrame = round(targetTime * fps);

    % Read frame
    video.CurrentTime = (targetFrame - 1) / fps;
    if hasFrame(video)
        frame = readFrame(video);

        % Calculate crop region (centered at detection position)
        x1 = max(1, cx - round(cropW/2));
        x2 = min(imgW, cx + round(cropW/2));
        y1 = max(1, cy - round(cropH/2));
        y2 = min(imgH, cy + round(cropH/2));

        % Crop image
        croppedFrame = frame(y1:y2, x1:x2, :);

        % Save image
        filename = sprintf('%s_t%ds.png', VIDEO_NAME, dt);
        filepath = fullfile(outputFolder, filename);
        imwrite(croppedFrame, filepath);

        fprintf('Saved: %s (frame %d, t=%.2fs)\n', filename, targetFrame, targetTime);
    else
        fprintf('Warning: Frame %d not available (t=%.2fs)\n', targetFrame, targetTime);
    end
end

fprintf('\nDone! Snapshots saved to: %s\n', outputFolder);
