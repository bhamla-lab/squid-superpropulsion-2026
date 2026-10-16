% C02: Batch process for green dye detection (Front view videos)
% 모든 Front view 동영상에 대해 C01을 개별 실행

clc; clear all; close all;

%% Define all Front view videos
% Format: {VIDEO_NAME, VIDEO_PATH}
videoList = {
    % Flex 5mm
    'F_F05mm_C01', '../01_rawdata/Flex5mm1mL2.5Vpump15Vactuator/FrontView/GX011411.MP4';
    'F_F05mm_C02', '../01_rawdata/Flex5mm1mL2.5Vpump15Vactuator/FrontView/GX011412.MP4';
    'F_F05mm_C03', '../01_rawdata/Flex5mm1mL2.5Vpump15Vactuator/FrontView/GX011413.MP4';
    % Flex 10mm
    'F_F10mm_C01', '../01_rawdata/Flex10mm1mL2.5Vpump15Vactuator/FrontView/GX011427.MP4';
    'F_F10mm_C02', '../01_rawdata/Flex10mm1mL2.5Vpump15Vactuator/FrontView/GX011428.MP4';
    'F_F10mm_C03', '../01_rawdata/Flex10mm1mL2.5Vpump15Vactuator/FrontView/GX011429.MP4';
    % Flex 15mm
    'F_F15mm_C01', '../01_rawdata/Flex15mm1mL2.5Vpump15Vactuator/FrontView/GX011416.MP4';
    'F_F15mm_C02', '../01_rawdata/Flex15mm1mL2.5Vpump15Vactuator/FrontView/GX011417.MP4';
    'F_F15mm_C03', '../01_rawdata/Flex15mm1mL2.5Vpump15Vactuator/FrontView/GX011421.MP4';
    % Flex 20mm
    'F_F20mm_C01', '../01_rawdata/Flex20mm1mL2.5Vpump15Vactuator/FrontView/GX011430.MP4';
    'F_F20mm_C02', '../01_rawdata/Flex20mm1mL2.5Vpump15Vactuator/FrontView/GX011431.MP4';
    'F_F20mm_C03', '../01_rawdata/Flex20mm1mL2.5Vpump15Vactuator/FrontView/GX011432.MP4';
    % Flex 30mm
    'F_F30mm_C01', '../01_rawdata/Flex30mm1mL2.5Vpump15Vactuator/FrontView/GX011422.MP4';
    'F_F30mm_C02', '../01_rawdata/Flex30mm1mL2.5Vpump15Vactuator/FrontView/GX011424.MP4';
    'F_F30mm_C03', '../01_rawdata/Flex30mm1mL2.5Vpump15Vactuator/FrontView/GX011425.MP4';
    'F_F30mm_C04', '../01_rawdata/Flex30mm1mL2.5Vpump15Vactuator/FrontView/GX011426.MP4';
    % Rigid 30mm
    'F_R30mm_C01', '../01_rawdata/Rigid1mL2.5Vpump15Vactuator/Front View/GX011408.MP4';
    'F_R30mm_C02', '../01_rawdata/Rigid1mL2.5Vpump15Vactuator/Front View/GX011409.MP4';
    'F_R30mm_C03', '../01_rawdata/Rigid1mL2.5Vpump15Vactuator/Front View/GX011410.MP4';
};

%% User selection
fprintf('=== Available Videos ===\n');
for i = 1:size(videoList, 1)
    matFile = fullfile('../03_output/00_start_detection', sprintf('%s.mat', videoList{i,1}));
    if exist(matFile, 'file')
        status = '[DONE]';
    else
        status = '[ ]';
    end
    fprintf('%d. %s %s\n', i, status, videoList{i,1});
end
fprintf('\n');

VIDEO_INDEX = input(sprintf('Select video index (1-%d): ', size(videoList, 1)));

if VIDEO_INDEX < 1 || VIDEO_INDEX > size(videoList, 1)
    error('Invalid index');
end

VIDEO_NAME = videoList{VIDEO_INDEX, 1};
VIDEO_PATH = videoList{VIDEO_INDEX, 2};

%% C01 Parameters
SAMPLE_INTERVAL = 1;
BRIGHTNESS_THRESHOLD = 50;  % Mean brightness < 50 = lights off
MOVING_STD_WINDOW = 5;      % Window size for moving std
STD_THRESHOLD = 1.0;        % Std threshold for green detection
window_size = 0.04;         % Detection window size (2% of frame)
EXTRA_FRAMES = 20;          % Continue monitoring after detection
MAT_OUTPUT_FOLDER = '../03_output/00_start_detection';

%% Check if already processed
matFileName = fullfile(MAT_OUTPUT_FOLDER, sprintf('%s.mat', VIDEO_NAME));
if exist(matFileName, 'file')
    overwrite = input(sprintf('%s already exists. Overwrite? (y/n): ', VIDEO_NAME), 's');
    if ~strcmpi(overwrite, 'y')
        fprintf('Cancelled.\n');
        return;
    end
end

%% Setup
if ~exist(MAT_OUTPUT_FOLDER, 'dir'), mkdir(MAT_OUTPUT_FOLDER); end

if ~exist(VIDEO_PATH, 'file')
    error('Video file not found: %s', VIDEO_PATH);
end

try
    video = VideoReader(VIDEO_PATH);
    fps = video.FrameRate;
catch ME
    fprintf('\n*** ERROR: Cannot read video file ***\n');
    fprintf('File: %s\n', VIDEO_PATH);
    fprintf('Error: %s\n', ME.message);
    fprintf('\nThis video may be corrupted or use an unsupported codec.\n');
    fprintf('Try re-encoding with: ffmpeg -i input.MP4 -c:v libx264 -c:a aac output.MP4\n');
    return;
end

fprintf('\n=== Processing: %s ===\n', VIDEO_NAME);
fprintf('Video: %s (%.2f fps)\n', video.Name, fps);
fprintf('Step 1: Waiting for lights to turn off...\n');

%% Step 1: Find when lights turn off
frameNum = 1;
darkStartFrame = [];

while hasFrame(video) && isempty(darkStartFrame)
    frame = readFrame(video);

    if mod(frameNum - 1, SAMPLE_INTERVAL) == 0
        meanBrightness = mean(rgb2gray(frame), 'all');

        if meanBrightness < BRIGHTNESS_THRESHOLD
            darkStartFrame = frameNum;
            fprintf('  -> Lights OFF at frame %d (%.2fs), brightness=%.1f\n', ...
                    frameNum, frameNum/fps, meanBrightness);
            break;
        end
    end
    frameNum = frameNum + 1;
end

if isempty(darkStartFrame)
    error('Lights never turned off in video');
end

%% Step 2: Find green dye after lights are off
fprintf('\nStep 2: Detecting green dye in center...\n');

greenStartFrame = [];
% Store green channel values over time
timeLog = [];
greenLog = [];
extraFrameCount = 0;

while hasFrame(video)
    frame = readFrame(video);

    if mod(frameNum - 1, SAMPLE_INTERVAL) == 0
        greenValue = detectGreenCenter(frame, window_size);

        % Log data
        timeLog(end+1) = frameNum / fps;
        greenLog(end+1) = greenValue;

        % Auto-detect using moving std
        if isempty(greenStartFrame) && length(greenLog) >= MOVING_STD_WINDOW
            recentStd = std(greenLog(end-MOVING_STD_WINDOW+1:end));
            if recentStd > STD_THRESHOLD
                greenStartFrame = frameNum;
                fprintf('  -> GREEN detected at frame %d (%.2fs), green=%.1f, std=%.2f\n', ...
                        frameNum, frameNum/fps, greenValue, recentStd);
            end
        end

        % Continue for EXTRA_FRAMES after detection
        if ~isempty(greenStartFrame)
            extraFrameCount = extraFrameCount + 1;
            if extraFrameCount >= EXTRA_FRAMES
                break;
            end
        end
    end
    frameNum = frameNum + 1;
end

%% Save results to mat file
% Set fixed detection center position
cx = 1884;
cy = 1079;

% Create data structure
startInfo.videoName = VIDEO_NAME;
startInfo.videoPath = VIDEO_PATH;
startInfo.fps = fps;
startInfo.darkFrame = darkStartFrame;
startInfo.darkTime = darkStartFrame / fps;
if ~isempty(greenStartFrame)
    startInfo.greenFrame = greenStartFrame;
    startInfo.greenTime = greenStartFrame / fps;
else
    startInfo.greenFrame = [];
    startInfo.greenTime = [];
end
startInfo.greenPixelPos = [cx, cy];

% Save to individual mat file in output folder
save(matFileName, 'startInfo');

%% Results
fprintf('\n=== Results ===\n');
fprintf('Lights OFF: frame %d (%.2f sec)\n', darkStartFrame, darkStartFrame/fps);
if ~isempty(greenStartFrame)
    fprintf('Green START: frame %d (%.2f sec)\n', greenStartFrame, greenStartFrame/fps);
    fprintf('Delay: %.2f sec\n', (greenStartFrame-darkStartFrame)/fps);
    fprintf('Detection center: [%d, %d]\n', cx, cy);
    fprintf('\nSaved to: %s\n', matFileName);
else
    fprintf('Green dye NOT detected\n');
end

%% Visualize before/after frames
video2 = VideoReader(VIDEO_PATH);

% Dark: before and after
darkBefore = max(1, darkStartFrame - SAMPLE_INTERVAL);
video2.CurrentTime = (darkBefore-1) / fps;
frameDarkBefore = readFrame(video2);

video2.CurrentTime = (darkStartFrame-1) / fps;
frameDarkAfter = readFrame(video2);

if ~isempty(greenStartFrame)
    % Green: before, detected, and after
    greenBefore = max(1, greenStartFrame - SAMPLE_INTERVAL);
    video2.CurrentTime = (greenBefore-1) / fps;
    frameGreenBefore = readFrame(video2);

    video2.CurrentTime = (greenStartFrame-1) / fps;
    frameGreenDetected = readFrame(video2);

    greenAfter = min(video.NumFrames, greenStartFrame + SAMPLE_INTERVAL);
    video2.CurrentTime = (greenAfter-1) / fps;
    frameGreenAfter = readFrame(video2);

    % Use fixed center position
    cx = 1884;
    cy = 1079;
    r = round(min(size(frameGreenDetected,2), size(frameGreenDetected,1)) * window_size);

    % Plot
    figure('Position', [100, 100, 2000, 1000]);

    subplot(2,3,1);
    imshow(frameDarkBefore);
    title(sprintf('Before Dark\nFrame %d (%.2fs)', darkBefore, darkBefore/fps), 'FontSize', 12);

    subplot(2,3,2);
    imshow(frameDarkAfter);
    title(sprintf('After Dark (Detected)\nFrame %d (%.2fs)', darkStartFrame, darkStartFrame/fps), 'FontSize', 12);

    subplot(2,3,3);
    imshow(frameGreenBefore); hold on;
    rectangle('Position', [cx-r, cy-r, 2*r, 2*r], 'EdgeColor', 'r', 'LineWidth', 2);
    title(sprintf('Before Green\nFrame %d (%.2fs)', greenBefore, greenBefore/fps), 'FontSize', 12);

    subplot(2,3,4);
    imshow(frameGreenDetected); hold on;
    rectangle('Position', [cx-r, cy-r, 2*r, 2*r], 'EdgeColor', 'r', 'LineWidth', 2);
    title(sprintf('Green Detected\nFrame %d (%.2fs)', greenStartFrame, greenStartFrame/fps), 'FontSize', 12);

    subplot(2,3,5);
    imshow(frameGreenAfter); hold on;
    rectangle('Position', [cx-r, cy-r, 2*r, 2*r], 'EdgeColor', 'r', 'LineWidth', 2);
    title(sprintf('After Green\nFrame %d (%.2fs)', greenAfter, greenAfter/fps), 'FontSize', 12);

    sgtitle(sprintf('%s | Delay: %.2f sec | Red box = Detection area', ...
            VIDEO_NAME, (greenStartFrame-darkStartFrame)/fps), 'FontSize', 14, 'FontWeight', 'bold');

    % Green channel monitoring plot
    figure('Position', [100, 100, 1200, 600]);
    plot(timeLog, greenLog, 'g-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    hold on;

    % Mark detection point
    detectionTime = greenStartFrame / fps;
    xline(detectionTime, 'r-', 'LineWidth', 3, 'Label', 'Auto-Detection');

    % Highlight detection region
    detectionIdx = find(timeLog == detectionTime);
    if ~isempty(detectionIdx)
        plot(timeLog(detectionIdx), greenLog(detectionIdx), 'ro', 'MarkerSize', 15, 'LineWidth', 3);
    end

    % Auto y-limits with 10% margin
    yMin = min(greenLog);
    yMax = max(greenLog);
    yMargin = (yMax - yMin) * 0.1;
    yLimits = [max(0, yMin - yMargin), min(255, yMax + yMargin)];

    % Shade the extra monitoring region (after detection)
    extraStartTime = detectionTime;
    extraEndTime = timeLog(end);
    fill([extraStartTime extraEndTime extraEndTime extraStartTime], ...
         [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
         'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Extra monitoring');

    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Green Channel Value (0-255)', 'FontSize', 12);
    title(sprintf('%s - Green Channel (Moving STD detection)', VIDEO_NAME), 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
    ylim(yLimits);
    grid off;
    legend('Location', 'best');
end

fprintf('\nDone! Close the figures to continue.\n');

%% Helper Function
function greenValue = detectGreenCenter(frame, window_size)
    % Fixed detection center
    cx = 1884;
    cy = 1079;
    [h, w, ~] = size(frame);
    r = round(min(w,h) * window_size);  % Center window size

    % Extract detection window (red box area only)
    centerFrame = frame(max(1,cy-r):min(h,cy+r), max(1,cx-r):min(w,cx+r), :);

    % Calculate mean green channel value (0-255)
    greenValue = mean(centerFrame(:,:,2), 'all');
end
