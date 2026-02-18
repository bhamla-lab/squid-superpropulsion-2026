% Binarize a specific frame from video
% === USER SETTINGS ===
videoFile = 'GX011456.MP4';   % Video filename (place in same folder or set full path)
frameNumber = 1500;      % Frame to extract and binarize
threshold = 0.5;        % Binarization threshold (0.0 to 1.0)
% ====================

videoPath = fullfile(fileparts(mfilename('fullpath')), videoFile);

% Read video
v = VideoReader(videoPath);

% Read specific frame
v.CurrentTime = (frameNumber - 1) / v.FrameRate;
frame = readFrame(v);

% Convert to grayscale
grayFrame = rgb2gray(frame);

% Load or create polygon ROI
roiFile = 'C02_polygon_roi.mat';
if exist(roiFile, 'file')
    load(roiFile, 'polyX', 'polyY');
    fprintf('Loaded existing polygon ROI\n');
else
    fprintf('Draw polygon ROI (double-click to finish)...\n');
    figure;
    imshow(grayFrame);
    title('Draw Polygon ROI - Double-click to finish');
    [polyX, polyY] = getline();
    close;
    save(roiFile, 'polyX', 'polyY');
    fprintf('Polygon ROI saved\n');
end

% Create mask
mask = poly2mask(polyX, polyY, size(grayFrame, 1), size(grayFrame, 2));

% Apply mask to grayscale
maskedGray = grayFrame;
maskedGray(~mask) = 0;

% Binarize and apply mask
binaryFrame = imbinarize(grayFrame, threshold);
binaryFrame(~mask) = 0;

% Display results
figure;
subplot(2, 2, 1);
imshow(frame);
hold on;
plot([polyX; polyX(1)], [polyY; polyY(1)], 'r-', 'LineWidth', 2);
hold off;
title(sprintf('Original Frame %d (ROI)', frameNumber));

subplot(2, 2, 2);
imshow(maskedGray);
title('Grayscale (ROI only)');

subplot(2, 2, 3);
imshow(binaryFrame);
title(sprintf('Binary (threshold=%.2f)', threshold));

subplot(2, 2, 4);
roiPixels = grayFrame(mask);
histogram(roiPixels, 256);
hold on;
yLim = ylim;
plot([threshold*255, threshold*255], yLim, 'r-', 'LineWidth', 2);
hold off;
title('Intensity Distribution (ROI)');
xlabel('Intensity');
ylabel('Pixel Count');
legend('Histogram', 'Threshold', 'Location', 'best');

% Save binary image
outputFolder = 'output';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
outputFile = fullfile(outputFolder, sprintf('binary_frame_%04d_th%.2f.png', frameNumber, threshold));
imwrite(binaryFrame, outputFile);

fprintf('Processed frame %d with threshold %.2f\n', frameNumber, threshold);
fprintf('Saved to: %s\n', outputFile);
