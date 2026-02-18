% Batch binarize frames from video
% === USER SETTINGS ===
videoFile = 'GX011456.MP4';   % Video filename (place in same folder or set full path)
frameInterval = 50;     % Process every Nth frame
threshold = 0.3;        % Binarization threshold (0.0 to 1.0)
% ====================

videoPath = fullfile(fileparts(mfilename('fullpath')), videoFile);

% Read video
v = VideoReader(videoPath);
totalFrames = floor(v.Duration * v.FrameRate);

% Load or create polygon ROI
roiFile = 'C02_polygon_roi.mat';
if exist(roiFile, 'file')
    load(roiFile, 'polyX', 'polyY');
    fprintf('Loaded existing polygon ROI\n');
else
    fprintf('Draw polygon ROI (double-click to finish)...\n');
    v.CurrentTime = 0;
    firstFrame = readFrame(v);
    grayFirst = rgb2gray(firstFrame);
    figure;
    imshow(grayFirst);
    title('Draw Polygon ROI - Double-click to finish');
    [polyX, polyY] = getline();
    close;
    save(roiFile, 'polyX', 'polyY');
    fprintf('Polygon ROI saved\n');
    v.CurrentTime = 0;
end

% Create output folder
outputFolder = 'output';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Process frames
frameCount = 0;
processedCount = 0;

fprintf('Processing frames with interval %d...\n', frameInterval);

while hasFrame(v)
    frame = readFrame(v);
    frameCount = frameCount + 1;

    if mod(frameCount, frameInterval) ~= 0
        continue;
    end

    % Convert to grayscale
    grayFrame = rgb2gray(frame);

    % Create mask
    mask = poly2mask(polyX, polyY, size(grayFrame, 1), size(grayFrame, 2));

    % Binarize and apply mask
    binaryFrame = imbinarize(grayFrame, threshold);
    binaryFrame(~mask) = 0;

    % Save
    outputFile = fullfile(outputFolder, sprintf('binary_frame_%04d_th%.2f.png', frameCount, threshold));
    imwrite(binaryFrame, outputFile);

    processedCount = processedCount + 1;

    % Progress
    progress = frameCount / totalFrames * 100;
    fprintf('\rProgress: %d/%d frames (%.1f%%) - %d saved', frameCount, totalFrames, progress, processedCount);
end
fprintf('\n');

fprintf('Batch processing complete: %d frames saved\n', processedCount);
