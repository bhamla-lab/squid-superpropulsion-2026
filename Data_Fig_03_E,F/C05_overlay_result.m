% Overlay C04 result with last frame from video
% === USER SETTINGS ===
clear all; clc; close all;

videoFile = 'GX011456.MP4';   % Video filename (place in same folder or set full path)
overlayAlpha = 0.9;  % Transparency (0.0 = only video, 1.0 = only overlay)
% ====================

videoPath = fullfile(fileparts(mfilename('fullpath')), videoFile);
overlayPath = fullfile('output', 'averaged_binary.png');     % C04 output

% Read last frame from video
v = VideoReader(videoPath);
v.CurrentTime = v.Duration - 0.1;
lastFrame = readFrame(v);

% Read overlay image (C04 result)
overlayImg = imread(overlayPath);

% Convert overlay to sky blue color (RGB: 135, 206, 235)
if size(overlayImg, 3) == 1
    skyBlueOverlay = zeros(size(overlayImg, 1), size(overlayImg, 2), 3, 'uint8');
    normalized = double(overlayImg) / 255;
    skyBlueOverlay(:, :, 1) = uint8(normalized * 135);  % Red
    skyBlueOverlay(:, :, 2) = uint8(normalized * 206);  % Green
    skyBlueOverlay(:, :, 3) = uint8(normalized * 235);  % Blue
    overlayImg = skyBlueOverlay;
end

% Ensure same size
if ~isequal(size(lastFrame), size(overlayImg))
    overlayImg = imresize(overlayImg, [size(lastFrame, 1), size(lastFrame, 2)]);
end

% Create overlay (raw image at full opacity, overlay with alpha)
blendedImg = double(lastFrame) + double(overlayImg) * overlayAlpha;
blendedImg(blendedImg > 255) = 255;  % Saturate
blendedImg = uint8(blendedImg);

% Display results in separate windows
figure('Name', 'Last Frame');
imshow(lastFrame);
title('Last Frame');

figure('Name', 'Overlay Result');
imshow(blendedImg);
title(sprintf('Overlay (alpha=%.2f)', overlayAlpha));

% Save result
outputFile = fullfile('output', 'overlay_result.png');
imwrite(blendedImg, outputFile);

fprintf('Overlay saved to: %s\n', outputFile);
