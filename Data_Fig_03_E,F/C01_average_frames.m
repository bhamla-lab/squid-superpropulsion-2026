% Average all frames from video
% === USER SETTINGS ===
videoFile = 'GX011456.MP4';   % Video filename (place in same folder or set full path)
skipFrames = 50;  % Process every Nth frame (1 = all frames, 10 = every 10th frame)
% ====================

videoPath = fullfile(fileparts(mfilename('fullpath')), videoFile);

% Read video
v = VideoReader(videoPath);
totalFrames = floor(v.Duration * v.FrameRate);

% Initialize accumulator
frameSum = [];
frameCount = 0;
skipCount = 0;

fprintf('Processing every %dth frame (%d total frames)...\n', skipFrames, totalFrames);

% Process frames
while hasFrame(v)
    frame = readFrame(v);

    skipCount = skipCount + 1;
    if mod(skipCount, skipFrames) ~= 0
        continue;
    end

    frame = double(frame);

    if isempty(frameSum)
        frameSum = frame;
    else
        frameSum = frameSum + frame;
    end

    frameCount = frameCount + 1;

    % Progress bar
    progress = skipCount / totalFrames * 100;
    fprintf('\rProgress: %d/%d frames (%.1f%%) - %d averaged', skipCount, totalFrames, progress, frameCount);
end
fprintf('\n');

% Calculate average
avgFrame = frameSum / frameCount;
avgFrame = uint8(avgFrame);

% Display result
figure;
imshow(avgFrame);
title(sprintf('Average of %d frames', frameCount));

% Save result
outputFolder = 'output';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
outputFile = fullfile(outputFolder, 'averaged_frames.png');
imwrite(avgFrame, outputFile);

fprintf('Averaged %d frames\n', frameCount);
fprintf('Saved to: %s\n', outputFile);
