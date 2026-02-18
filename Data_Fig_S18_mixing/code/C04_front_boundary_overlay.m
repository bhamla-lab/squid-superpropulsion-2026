% C04: Front view - Overlay boundaries for Rigid vs Flexible

clc; clear all; close all;

%% Parameters
SNAPSHOT_BASE = '../data/snapshots';
NOZZLE_DIAMETER_PX = 13.5;  % Nozzle diameter in pixels

rigidVideo = 'F_R30mm_C01';
flexVideo = 'F_F20mm_C01';

fprintf('Comparing: %s vs %s\n', rigidVideo, flexVideo);

%% Extract boundaries from first frames
rigidFolder = fullfile(SNAPSHOT_BASE, rigidVideo);
flexFolder = fullfile(SNAPSHOT_BASE, flexVideo);

rigidImages = dir(fullfile(rigidFolder, '*.png'));
flexImages = dir(fullfile(flexFolder, '*.png'));

[~, sortIdx] = sort(cellfun(@(x) extractTimeFromFilename(x), {rigidImages.name}));
rigidImages = rigidImages(sortIdx);
[~, sortIdx] = sort(cellfun(@(x) extractTimeFromFilename(x), {flexImages.name}));
flexImages = flexImages(sortIdx);

numFrames = min(4, min(length(rigidImages), length(flexImages)));

rigidBoundaries = {};
for i = 1:numFrames
    img = imread(fullfile(rigidFolder, rigidImages(i).name));
    [~, boundary, ~, ~] = detectGreenBoundary(img);
    if ~isempty(boundary)
        rigidBoundaries{end+1} = boundary;
    end
end

flexBoundaries = {};
for i = 1:numFrames
    img = imread(fullfile(flexFolder, flexImages(i).name));
    [~, boundary, ~, ~] = detectGreenBoundary(img);
    if ~isempty(boundary)
        flexBoundaries{end+1} = boundary;
    end
end

%% Plot overlaid boundaries
img = imread(fullfile(rigidFolder, rigidImages(1).name));
[imgH, imgW, ~] = size(img);
cx = imgW / 2;
cy = imgH / 2;

% Rigid (Blue)
fig1 = figure('Position', [100, 100, 600, 500]);
hold on;
blueColors = [linspace(0, 0, length(rigidBoundaries))', ...
              linspace(0, 0, length(rigidBoundaries))', ...
              linspace(0.3, 1, length(rigidBoundaries))'];
for i = 1:length(rigidBoundaries)
    x_norm = (rigidBoundaries{i}(:,2) - cx) / NOZZLE_DIAMETER_PX;
    y_norm = (rigidBoundaries{i}(:,1) - cy) / NOZZLE_DIAMETER_PX;
    fill(x_norm, y_norm, blueColors(i,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
text(-35, 25, sprintf('t=0~%ds', (numFrames-1)*2), 'FontSize', 10, 'BackgroundColor', 'w');
xlabel('x / D_{nozzle}', 'FontSize', 11);
ylabel('y / D_{nozzle}', 'FontSize', 11);
title(sprintf('Rigid (%s)', rigidVideo), 'FontSize', 12, 'Interpreter', 'none');
xlim([-40 50]); ylim([-30 30]);
axis equal; grid off;

% Flexible (Red)
fig2 = figure('Position', [100, 100, 600, 500]);
hold on;
redColors = [linspace(0.3, 1, length(flexBoundaries))', ...
             linspace(0, 0, length(flexBoundaries))', ...
             linspace(0, 0, length(flexBoundaries))'];
for i = 1:length(flexBoundaries)
    x_norm = (flexBoundaries{i}(:,2) - cx) / NOZZLE_DIAMETER_PX;
    y_norm = (flexBoundaries{i}(:,1) - cy) / NOZZLE_DIAMETER_PX;
    fill(x_norm, y_norm, redColors(i,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
text(-35, 25, sprintf('t=0~%ds', (numFrames-1)*2), 'FontSize', 10, 'BackgroundColor', 'w');
xlabel('x / D_{nozzle}', 'FontSize', 11);
ylabel('y / D_{nozzle}', 'FontSize', 11);
title(sprintf('Flexible (%s)', flexVideo), 'FontSize', 12, 'Interpreter', 'none');
xlim([-40 50]); ylim([-30 30]);
axis equal; grid off;

fprintf('\nVisualization complete (Front View)!\n');

%% Helper Functions
function timeVal = extractTimeFromFilename(filename)
    tokens = regexp(filename, '_t(\d+)s', 'tokens');
    if ~isempty(tokens)
        timeVal = str2double(tokens{1}{1});
    else
        timeVal = 0;
    end
end

function [greenMask, boundary, maxArea, avgIntensity] = detectGreenBoundary(img)
    hsv = rgb2hsv(img);
    greenMask = (hsv(:,:,1) >= 0.2 & hsv(:,:,1) <= 0.5) & ...
                (hsv(:,:,2) >= 0.2) & (hsv(:,:,3) >= 0.2);
    if any(greenMask(:))
        greenMask = bwareaopen(greenMask, 50);
        greenMask = imfill(greenMask, 'holes');
        boundaries = bwboundaries(greenMask, 'noholes');
        if ~isempty(boundaries)
            maxArea = 0; maxIdx = 1;
            for i = 1:length(boundaries)
                area = polyarea(boundaries{i}(:,2), boundaries{i}(:,1));
                if area > maxArea
                    maxArea = area; maxIdx = i;
                end
            end
            boundary = boundaries{maxIdx};
            avgIntensity = mean(hsv(:,:,2), 'all');
        else
            boundary = []; maxArea = 0; avgIntensity = 0;
        end
    else
        boundary = []; maxArea = 0; avgIntensity = 0;
    end
end
