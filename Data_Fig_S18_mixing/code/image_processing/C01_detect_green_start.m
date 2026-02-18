% C01: Detect when green dye appears after lights turn off
% 불이 꺼지고 암전된 후 중앙에 초록색이 나타나는 순간 찾기
clc; clear all; close all;

%% Parameters
VIDEO_NAME = 'S_F05mm_C01';  % Output mat file name
VIDEO_PATH = '../01_rawdata/Flex5mm1mL2.5Vpump15Vactuator/FrontView/GX011411.MP4';
SAMPLE_INTERVAL = 1;
BRIGHTNESS_THRESHOLD = 50;  % Mean brightness < 50 = lights off
GREEN_THRESHOLD = 0.01;      % 0.01% green in center = dye detected
MAT_OUTPUT_FOLDER = '../03_output/00_start_detection';  % Folder to save mat file

%% Setup
if ~exist(MAT_OUTPUT_FOLDER, 'dir'), mkdir(MAT_OUTPUT_FOLDER); end

video = VideoReader(VIDEO_PATH);
fps = video.FrameRate;

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
while hasFrame(video)
    frame = readFrame(video);

    if mod(frameNum - 1, SAMPLE_INTERVAL) == 0
        greenPercent = detectGreenCenter(frame);

        if greenPercent > GREEN_THRESHOLD
            greenStartFrame = frameNum;
            fprintf('  -> GREEN detected at frame %d (%.2fs), green=%.2f%%\n', ...
                    frameNum, frameNum/fps, greenPercent);
            break;
        end
    end
    frameNum = frameNum + 1;
end

%% Save results to mat file
% Get image dimensions and center position
tempVideo = VideoReader(VIDEO_PATH);
tempFrame = read(tempVideo, 1);
[h, w, ~] = size(tempFrame);
cx = round(w/2);
cy = round(h/2);

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
matFileName = fullfile(MAT_OUTPUT_FOLDER, sprintf('%s.mat', VIDEO_NAME));
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
    % Green: before and after
    greenBefore = max(1, greenStartFrame - SAMPLE_INTERVAL);
    video2.CurrentTime = (greenBefore-1) / fps;
    frameGreenBefore = readFrame(video2);

    video2.CurrentTime = (greenStartFrame-1) / fps;
    frameGreenAfter = readFrame(video2);

    % Calculate center region box (30% of image)
    [h, w, ~] = size(frameGreenAfter);
    cx = round(w/2); cy = round(h/2);
    r = round(min(w,h) * 0.15);  % 30% diameter = 15% radius

    % Plot
    figure('Position', [100, 100, 1600, 800]);

    subplot(2,2,1);
    imshow(frameDarkBefore);
    title(sprintf('Before Dark\nFrame %d (%.2fs)', darkBefore, darkBefore/fps), 'FontSize', 12);

    subplot(2,2,2);
    imshow(frameDarkAfter);
    title(sprintf('After Dark (Detected)\nFrame %d (%.2fs)', darkStartFrame, darkStartFrame/fps), 'FontSize', 12);

    subplot(2,2,3);
    imshow(frameGreenBefore); hold on;
    rectangle('Position', [cx-r, cy-r, 2*r, 2*r], 'EdgeColor', 'r', 'LineWidth', 2);
    title(sprintf('Before Green\nFrame %d (%.2fs)', greenBefore, greenBefore/fps), 'FontSize', 12);

    subplot(2,2,4);
    imshow(frameGreenAfter); hold on;
    rectangle('Position', [cx-r, cy-r, 2*r, 2*r], 'EdgeColor', 'r', 'LineWidth', 2);
    title(sprintf('After Green (Detected)\nFrame %d (%.2fs)', greenStartFrame, greenStartFrame/fps), 'FontSize', 12);

    sgtitle(sprintf('Detection Results | Delay: %.2f sec | Red box = Detection area (center 30%%)', ...
            (greenStartFrame-darkStartFrame)/fps), 'FontSize', 14, 'FontWeight', 'bold');
end

%% Helper Functions
function greenPercent = detectGreenCenter(frame)
    [h, w, ~] = size(frame);
    cx = round(w/2); cy = round(h/2);
    r = round(min(w,h) * 0.15);  % Center 30% diameter

    centerFrame = frame(max(1,cy-r):min(h,cy+r), max(1,cx-r):min(w,cx+r), :);
    hsv = rgb2hsv(centerFrame);

    greenMask = (hsv(:,:,1) >= 0.2 & hsv(:,:,1) <= 0.5) & ...
                (hsv(:,:,2) >= 0.2) & (hsv(:,:,3) >= 0.2);

    greenPercent = sum(greenMask(:)) / numel(greenMask) * 100;
end
