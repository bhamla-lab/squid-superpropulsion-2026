% C03: Compare normalized area at tau/T=400 vs nozzle length

clc; clear all; close all;

%% Parameters
METRICS_BASE = '../data/metrics';
NOZZLE_DIAMETER_PX = 13.5;  % Nozzle diameter in pixels
JET_VELOCITY_PX_SEC = 616.22;  % Jet velocity in px/sec
TARGET_TIME = 400;  % Target normalized time tau/T

%% Nozzle length to tau/T mapping
nozzleLengthMap = {
    'F05mm', 0.076980036;
    'F10mm', 0.153960072;
    'F15mm', 0.230940108;
    'F20mm', 0.307920144;
    'F30mm', 0.461880215;
    'R30mm', 0;
};

%% Get all metrics mat files
allDirs = dir(METRICS_BASE);
allDirs = allDirs([allDirs.isdir] & ~ismember({allDirs.name}, {'.', '..'}));

if isempty(allDirs)
    error('No metrics files found.');
end

fprintf('Found %d metrics files\n\n', length(allDirs));

%% Load and organize data by nozzle length
dataByNozzle = struct();  % Front view data
dataByNozzle_side = struct();  % Side view data

for i = 1:length(allDirs)
    matFile = fullfile(METRICS_BASE, allDirs(i).name, sprintf('%s_metrics.mat', allDirs(i).name));
    if exist(matFile, 'file')
        load(matFile, 'timeResults', 'VIDEO_NAME');

        isSideView = startsWith(VIDEO_NAME, 'S_');

        tokens = regexp(VIDEO_NAME, '_([FR]\d+mm)_C(\d+)', 'tokens');
        if ~isempty(tokens)
            nozzleKey = tokens{1}{1};
            caseNum = str2double(tokens{1}{2});

            % Normalize time: tau/T = t * U_jet / D_nozzle
            times = [timeResults.time] * JET_VELOCITY_PX_SEC / NOZZLE_DIAMETER_PX;

            nozzleArea = pi * (NOZZLE_DIAMETER_PX/2)^2;
            areas = [timeResults.area] / nozzleArea;
            intensities = [timeResults.intensity];

            % Find area at tau/T=400
            if ~isempty(times) && max(times) >= TARGET_TIME
                areaAtTarget = interp1(times, areas, TARGET_TIME, 'linear', NaN);

                if ~isnan(areaAtTarget)
                    if isSideView
                        if ~isfield(dataByNozzle_side, nozzleKey)
                            dataByNozzle_side.(nozzleKey).areas = [];
                            dataByNozzle_side.(nozzleKey).intensities = [];
                        end
                        dataByNozzle_side.(nozzleKey).areas(end+1) = areaAtTarget;
                    else
                        if ~isfield(dataByNozzle, nozzleKey)
                            dataByNozzle.(nozzleKey).areas = [];
                            dataByNozzle.(nozzleKey).intensities = [];
                        end
                        dataByNozzle.(nozzleKey).areas(end+1) = areaAtTarget;
                    end
                    fprintf('Loaded: %s - Area at tau/T=%d: %.2f\n', VIDEO_NAME, TARGET_TIME, areaAtTarget);
                else
                    fprintf('Warning: %s - tau/T=%d outside range\n', VIDEO_NAME, TARGET_TIME);
                end
            else
                fprintf('Warning: %s - insufficient data for tau/T=%d\n', VIDEO_NAME, TARGET_TIME);
            end

            % Find intensity at last time frame
            if ~isempty(times)
                lastIntensity = intensities(end);
                if isSideView
                    if ~isfield(dataByNozzle_side, nozzleKey)
                        dataByNozzle_side.(nozzleKey).areas = [];
                        dataByNozzle_side.(nozzleKey).intensities = [];
                    end
                    dataByNozzle_side.(nozzleKey).intensities(end+1) = lastIntensity;
                else
                    if ~isfield(dataByNozzle, nozzleKey)
                        dataByNozzle.(nozzleKey).areas = [];
                        dataByNozzle.(nozzleKey).intensities = [];
                    end
                    dataByNozzle.(nozzleKey).intensities(end+1) = lastIntensity;
                end
                fprintf('         %s - Intensity at last frame (tau/T=%.1f): %.4f\n', VIDEO_NAME, times(end), lastIntensity);
            end
        end
    end
end

%% Calculate statistics and prepare plot data
% Front view
nozzleLengths = [];
meanAreas = [];
stdAreas = [];
meanIntensities = [];
stdIntensities = [];
labels = {};

for i = 1:size(nozzleLengthMap, 1)
    nozzleKey = nozzleLengthMap{i, 1};
    nozzleLength_tauT = nozzleLengthMap{i, 2};

    if isfield(dataByNozzle, nozzleKey) && ~isempty(dataByNozzle.(nozzleKey).areas)
        nozzleLengths(end+1) = nozzleLength_tauT;
        meanAreas(end+1) = mean(dataByNozzle.(nozzleKey).areas);
        stdAreas(end+1) = std(dataByNozzle.(nozzleKey).areas);
        meanIntensities(end+1) = mean(dataByNozzle.(nozzleKey).intensities);
        stdIntensities(end+1) = std(dataByNozzle.(nozzleKey).intensities);
        labels{end+1} = nozzleKey;

        fprintf('\nFront View - %s (tau/T=%.3f):\n', nozzleKey, nozzleLength_tauT);
        fprintf('  Area: Mean=%.2f, Std=%.2f (n=%d)\n', ...
                meanAreas(end), stdAreas(end), length(dataByNozzle.(nozzleKey).areas));
        fprintf('  Intensity: Mean=%.4f, Std=%.4f (n=%d)\n', ...
                meanIntensities(end), stdIntensities(end), length(dataByNozzle.(nozzleKey).intensities));
    end
end

% Side view
nozzleLengths_side = [];
meanAreas_side = [];
stdAreas_side = [];
meanIntensities_side = [];
stdIntensities_side = [];
labels_side = {};

for i = 1:size(nozzleLengthMap, 1)
    nozzleKey = nozzleLengthMap{i, 1};
    nozzleLength_tauT = nozzleLengthMap{i, 2};

    if isfield(dataByNozzle_side, nozzleKey) && ~isempty(dataByNozzle_side.(nozzleKey).areas)
        nozzleLengths_side(end+1) = nozzleLength_tauT;
        meanAreas_side(end+1) = mean(dataByNozzle_side.(nozzleKey).areas);
        stdAreas_side(end+1) = std(dataByNozzle_side.(nozzleKey).areas);
        meanIntensities_side(end+1) = mean(dataByNozzle_side.(nozzleKey).intensities);
        stdIntensities_side(end+1) = std(dataByNozzle_side.(nozzleKey).intensities);
        labels_side{end+1} = nozzleKey;

        fprintf('\nSide View - %s (tau/T=%.3f):\n', nozzleKey, nozzleLength_tauT);
        fprintf('  Area: Mean=%.2f, Std=%.2f (n=%d)\n', ...
                meanAreas_side(end), stdAreas_side(end), length(dataByNozzle_side.(nozzleKey).areas));
        fprintf('  Intensity: Mean=%.4f, Std=%.4f (n=%d)\n', ...
                meanIntensities_side(end), stdIntensities_side(end), length(dataByNozzle_side.(nozzleKey).intensities));
    end
end

%% Plot
figure('Position', [100, 100, 2000, 600]);

% Subplot 1: Area at tau/T=400 (Front View)
subplot(1,3,1);
errorbar(nozzleLengths, meanAreas, stdAreas, 'o', 'LineWidth', 2, 'MarkerSize', 10, ...
         'CapSize', 10);
hold on;

for i = 1:length(labels)
    text(nozzleLengths(i), meanAreas(i) + stdAreas(i) + 5, labels{i}, ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end

xlabel('Nozzle Length (\tau_{nozzle}/T)', 'FontSize', 12);
ylabel(sprintf('Normalized Area at \\tau/T=%d (A/A_{nozzle})', TARGET_TIME), 'FontSize', 12);
title(sprintf('Front View: Area at \\tau/T=%d', TARGET_TIME), 'FontSize', 14, 'FontWeight', 'bold');
xlim([-0.05 0.5]);
grid off;

% Subplot 2: Area at tau/T=400 (Side View)
subplot(1,3,2);
errorbar(nozzleLengths_side, meanAreas_side, stdAreas_side, 'o', 'LineWidth', 2, 'MarkerSize', 10, ...
         'CapSize', 10, 'Color', [0.85 0.33 0.10]);
hold on;

for i = 1:length(labels_side)
    text(nozzleLengths_side(i), meanAreas_side(i) + stdAreas_side(i) + 5, labels_side{i}, ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end

xlabel('Nozzle Length (\tau_{nozzle}/T)', 'FontSize', 12);
ylabel(sprintf('Normalized Area at \\tau/T=%d (A/A_{nozzle})', TARGET_TIME), 'FontSize', 12);
title(sprintf('Side View: Area at \\tau/T=%d', TARGET_TIME), 'FontSize', 14, 'FontWeight', 'bold');
xlim([-0.05 0.5]);
grid off;

% Subplot 3: Intensity at last time frame
subplot(1,3,3);
errorbar(nozzleLengths, meanIntensities, stdIntensities, 's', 'LineWidth', 2, 'MarkerSize', 10, ...
         'MarkerFaceColor', 'r', 'CapSize', 10);
hold on;
xlim([-0.05 0.5]);

for i = 1:length(labels)
    text(nozzleLengths(i), meanIntensities(i) + stdIntensities(i) + 0.005, labels{i}, ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end

xlabel('Nozzle Length (\tau_{nozzle}/T)', 'FontSize', 12);
ylabel('Intensity at Last Frame (Saturation)', 'FontSize', 12);
title('Front View: Final Intensity', 'FontSize', 14, 'FontWeight', 'bold');
grid off;

sgtitle('Nozzle Length Effects on Mixing', 'FontSize', 16, 'FontWeight', 'bold');

%% Bar chart: Rigid vs Best Flexible
figure('Position', [100, 100, 800, 600]);

rigidIdx = find(strcmp(labels, 'R30mm'));
rigidArea = meanAreas(rigidIdx);
rigidAreaStd = stdAreas(rigidIdx);
rigidIntensity = meanIntensities(rigidIdx);
rigidIntensityStd = stdIntensities(rigidIdx);

flexIdx = find(~strcmp(labels, 'R30mm'));
[bestFlexArea, bestIdx] = max(meanAreas(flexIdx));
bestFlexAreaStd = stdAreas(flexIdx(bestIdx));
bestFlexIntensity = meanIntensities(flexIdx(bestIdx));
bestFlexIntensityStd = stdIntensities(flexIdx(bestIdx));
bestFlexLabel = labels{flexIdx(bestIdx)};

subplot(1,2,1);
b = bar([1 2], [rigidArea, bestFlexArea]);
b.FaceColor = 'flat';
b.CData(1,:) = [0.8 0.2 0.2];
b.CData(2,:) = [0.2 0.6 0.8];
hold on;
errorbar([1 2], [rigidArea, bestFlexArea], [rigidAreaStd, bestFlexAreaStd], 'k.', 'LineWidth', 1.5, 'CapSize', 10);
set(gca, 'XTickLabel', {'Rigid (R30mm)', sprintf('Flexible (%s)', bestFlexLabel)});
ylabel(sprintf('Normalized Area at \\tau/T=%d', TARGET_TIME), 'FontSize', 12);

areaIncrease = ((bestFlexArea - rigidArea) / rigidArea) * 100;
if areaIncrease > 0
    text(1.5, max([rigidArea, bestFlexArea]) * 1.1, sprintf('+%.1f%%', areaIncrease), ...
         'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0 0.6 0]);
else
    text(1.5, max([rigidArea, bestFlexArea]) * 1.1, sprintf('%.1f%%', areaIncrease), ...
         'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r');
end

title('Area Comparison', 'FontSize', 14, 'FontWeight', 'bold');
grid off;

subplot(1,2,2);
b = bar([1 2], [rigidIntensity, bestFlexIntensity]);
b.FaceColor = 'flat';
b.CData(1,:) = [0.8 0.2 0.2];
b.CData(2,:) = [0.2 0.6 0.8];
hold on;
errorbar([1 2], [rigidIntensity, bestFlexIntensity], [rigidIntensityStd, bestFlexIntensityStd], 'k.', 'LineWidth', 1.5, 'CapSize', 10);
set(gca, 'XTickLabel', {'Rigid (R30mm)', sprintf('Flexible (%s)', bestFlexLabel)});
ylabel('Intensity at Last Frame', 'FontSize', 12);

intensityIncrease = ((bestFlexIntensity - rigidIntensity) / rigidIntensity) * 100;
if intensityIncrease > 0
    text(1.5, max([rigidIntensity, bestFlexIntensity]) * 1.1, sprintf('+%.1f%%', intensityIncrease), ...
         'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0 0.6 0]);
else
    text(1.5, max([rigidIntensity, bestFlexIntensity]) * 1.1, sprintf('%.1f%%', intensityIncrease), ...
         'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r');
end

title('Intensity Comparison', 'FontSize', 14, 'FontWeight', 'bold');
grid off;

sgtitle('Rigid vs Best Flexible Nozzle', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('\n=== Comparison Results ===\n');
fprintf('Best flexible nozzle: %s (Area=%.2f+/-%.2f, Intensity=%.4f+/-%.4f)\n', ...
        bestFlexLabel, bestFlexArea, bestFlexAreaStd, bestFlexIntensity, bestFlexIntensityStd);
fprintf('Rigid nozzle: R30mm (Area=%.2f+/-%.2f, Intensity=%.4f+/-%.4f)\n', ...
        rigidArea, rigidAreaStd, rigidIntensity, rigidIntensityStd);
fprintf('\nArea increase: %.1f%%\n', areaIncrease);
fprintf('Intensity increase: %.1f%%\n', intensityIncrease);

fprintf('\nVisualization complete!\n');
