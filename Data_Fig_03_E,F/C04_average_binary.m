% Average all binary images in output folder
% === USER SETTINGS ===
clear all; clc; close all;

weightFactor = 10;  % Weight for white pixels (higher = more emphasis on white)
lowLimit = 1;      % Lower intensity limit (0-255)
upperLimit = 20;   % Upper intensity limit (0-255)
% ====================

outputFolder = 'output';

% Get all PNG files
imageFiles = dir(fullfile(outputFolder, '*.png'));
numImages = length(imageFiles);

if numImages == 0
    error('No images found in output folder');
end

fprintf('Averaging %d images...\n', numImages);

% Initialize accumulator
imageSum = [];

% Process all images
for i = 1:numImages
    imgPath = fullfile(outputFolder, imageFiles(i).name);
    img = imread(imgPath);
    img = double(img) * weightFactor;  % Weight white pixels

    if isempty(imageSum)
        imageSum = img;
    else
        imageSum = imageSum + img;
    end

    fprintf('\rProgress: %d/%d images (%.1f%%)', i, numImages, i/numImages*100);
end
fprintf('\n');

% Calculate weighted average
avgImage = imageSum / numImages;
avgImage = avgImage / max(avgImage(:)) * 255;  % Scale to 0-255
avgImage = uint8(avgImage);

% Apply intensity clipping and rescale to 0-255
clippedImage = double(avgImage);
clippedImage(clippedImage < lowLimit) = lowLimit;
clippedImage(clippedImage > upperLimit) = upperLimit;
clippedImage = (clippedImage - lowLimit) / (upperLimit - lowLimit) * 255;
clippedImage = uint8(clippedImage);

% Display results
figure;
subplot(2, 2, 1);
imshow(avgImage);
title('Original Average');

subplot(2, 2, 2);
imshow(clippedImage);
title(sprintf('Clipped [%d, %d]', lowLimit, upperLimit));

subplot(2, 2, 3);
histogram(avgImage, 256);
hold on;
yLim = ylim;
plot([lowLimit, lowLimit], yLim, 'r-', 'LineWidth', 2);
plot([upperLimit, upperLimit], yLim, 'g-', 'LineWidth', 2);
hold off;
title('Intensity Distribution');
xlabel('Intensity');
ylabel('Pixel Count');
legend('Histogram', 'Low Limit', 'Upper Limit', 'Location', 'best');

subplot(2, 2, 4);
histogram(clippedImage, 256);
title('Clipped Distribution');
xlabel('Intensity');
ylabel('Pixel Count');

% Save result
outputFile = fullfile(outputFolder, 'averaged_binary.png');
imwrite(clippedImage, outputFile);

fprintf('Averaged %d images\n', numImages);
fprintf('Saved to: %s\n', outputFile);
