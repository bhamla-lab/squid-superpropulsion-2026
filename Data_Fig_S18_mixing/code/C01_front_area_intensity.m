% C01: Front view - Area and Intensity vs Normalized Time

clc; clear all; close all;

%% Parameters
METRICS_BASE = '../data/metrics';
NOZZLE_DIAMETER_PX = 13.5;  % Nozzle diameter in pixels
JET_VELOCITY_PX_SEC = 616.22;  % Jet velocity in px/sec

%% Nozzle length mapping
nozzleLengthMap = containers.Map(...
    {'F05mm', 'F10mm', 'F15mm', 'F20mm', 'F30mm', 'R30mm'}, ...
    {0.076980036, 0.153960072, 0.230940108, 0.307920144, 0.461880215, 0});

nozzleDisplayNames = containers.Map(...
    {'F05mm', 'F10mm', 'F15mm', 'F20mm', 'F30mm', 'R30mm'}, ...
    {'Flex (\tau/T=0.077)', 'Flex (\tau/T=0.154)', 'Flex (\tau/T=0.231)', 'Flex (\tau/T=0.308)', 'Flex (\tau/T=0.462)', 'Rigid (\tau/T=0)'});

%% Get all metrics mat files
allDirs = dir(METRICS_BASE);
allDirs = allDirs([allDirs.isdir] & ~ismember({allDirs.name}, {'.', '..'}));

allMetrics = {};
for i = 1:length(allDirs)
    matFile = fullfile(METRICS_BASE, allDirs(i).name, sprintf('%s_metrics.mat', allDirs(i).name));
    if exist(matFile, 'file')
        allMetrics{end+1} = matFile;
    end
end

if isempty(allMetrics)
    error('No metrics files found.');
end

fprintf('Found %d metrics files\n\n', length(allMetrics));

%% Load and organize data by chamber size
dataBySize = struct();

for i = 1:length(allMetrics)
    load(allMetrics{i}, 'timeResults', 'VIDEO_NAME');

    tokens = regexp(VIDEO_NAME, '_([FR])(\d+)mm_C(\d+)', 'tokens');
    if ~isempty(tokens)
        chamberType = tokens{1}{1};
        chamberSize = str2double(tokens{1}{2});

        sizeKey = sprintf('%s%02dmm', chamberType, chamberSize);

        if ~isfield(dataBySize, sizeKey)
            dataBySize.(sizeKey) = struct('videos', {}, 'data', {});
        end

        idx = length(dataBySize.(sizeKey)) + 1;
        dataBySize.(sizeKey)(idx).videoName = VIDEO_NAME;
        dataBySize.(sizeKey)(idx).times = [timeResults.time] * JET_VELOCITY_PX_SEC / NOZZLE_DIAMETER_PX;
        nozzleArea = pi * (NOZZLE_DIAMETER_PX/2)^2;
        dataBySize.(sizeKey)(idx).areas = [timeResults.area] / nozzleArea;
        dataBySize.(sizeKey)(idx).intensities = [timeResults.intensity];

        if isKey(nozzleLengthMap, sizeKey)
            dataBySize.(sizeKey)(idx).nozzleLength = nozzleLengthMap(sizeKey);
            dataBySize.(sizeKey)(idx).displayName = nozzleDisplayNames(sizeKey);
        else
            dataBySize.(sizeKey)(idx).nozzleLength = NaN;
            dataBySize.(sizeKey)(idx).displayName = sizeKey;
        end

        [sortedTimes, sortIdx] = sort(dataBySize.(sizeKey)(idx).times);
        dataBySize.(sizeKey)(idx).times = sortedTimes;
        dataBySize.(sizeKey)(idx).areas = dataBySize.(sizeKey)(idx).areas(sortIdx);
        dataBySize.(sizeKey)(idx).intensities = dataBySize.(sizeKey)(idx).intensities(sortIdx);

        fprintf('Loaded: %s (%s, tau/T=%.3f)\n', VIDEO_NAME, sizeKey, dataBySize.(sizeKey)(idx).nozzleLength);
    end
end

%% Plot: Area and Intensity vs Time
figure('Position', [100, 100, 1400, 600]);

sizeNames = fieldnames(dataBySize);
colors = lines(length(sizeNames));

for i = 1:length(sizeNames)
    sizeName = sizeNames{i};
    data = dataBySize.(sizeName);

    allTimes = [];
    allAreas = [];

    for j = 1:length(data)
        allTimes = [allTimes; data(j).times(:)];
        allAreas = [allAreas; data(j).areas(:)];
    end

    subplot(1,2,1); hold on;
    for j = 1:length(data)
        plot(data(j).times, data(j).areas, 'o-', 'Color', colors(i,:), 'LineWidth', 1, ...
             'MarkerSize', 4, 'HandleVisibility', 'off');
    end

    uniqueTimes = unique(allTimes);
    meanAreas = arrayfun(@(t) mean(allAreas(allTimes == t)), uniqueTimes);
    displayName = data(1).displayName;
    plot(uniqueTimes, meanAreas, 'o-', 'Color', colors(i,:), 'LineWidth', 2, ...
         'MarkerSize', 8, 'DisplayName', displayName, 'MarkerFaceColor', colors(i,:));
end

subplot(1,2,1);
xlabel('Normalized Time (t* = t \cdot U_{jet} / D_{nozzle})', 'FontSize', 12);
ylabel('Normalized Area (A/A_{nozzle})', 'FontSize', 12);
title('Green Dye Area vs Normalized Time', 'FontSize', 14);
legend('Location', 'best');
xlim([0 1000]);
grid off;

for i = 1:length(sizeNames)
    sizeName = sizeNames{i};
    data = dataBySize.(sizeName);

    allTimes = [];
    allIntensities = [];

    for j = 1:length(data)
        allTimes = [allTimes; data(j).times(:)];
        allIntensities = [allIntensities; data(j).intensities(:)];
    end

    subplot(1,2,2); hold on;
    for j = 1:length(data)
        plot(data(j).times, data(j).intensities, 'o-', 'Color', colors(i,:), 'LineWidth', 1, ...
             'MarkerSize', 4, 'HandleVisibility', 'off');
    end

    uniqueTimes = unique(allTimes);
    meanIntensities = arrayfun(@(t) mean(allIntensities(allTimes == t)), uniqueTimes);
    displayName = data(1).displayName;
    plot(uniqueTimes, meanIntensities, 'o-', 'Color', colors(i,:), 'LineWidth', 2, ...
         'MarkerSize', 8, 'DisplayName', displayName, 'MarkerFaceColor', colors(i,:));
end

subplot(1,2,2);
xlabel('Normalized Time (t* = t \cdot U_{jet} / D_{nozzle})', 'FontSize', 12);
ylabel('Intensity (Saturation)', 'FontSize', 12);
title('Green Dye Intensity vs Normalized Time', 'FontSize', 14);
legend('Location', 'best');
grid on;

sgtitle('Front View: Mixing Metrics by Nozzle Type', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('\nVisualization complete!\n');
