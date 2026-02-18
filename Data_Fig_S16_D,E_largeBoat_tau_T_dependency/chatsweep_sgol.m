% ChatSweep Data Analysis using Savitzky-Golay smoothing
% Compares different nozzle lengths at 9V

% METHOD:
%   1. Applies Savitzky-Golay filter for smoothing (preserves peaks/features)
%   2. Calculates velocity using numerical differentiation (gradient)
%   3. Calculates acceleration by differentiating velocity
%   4. Reports smoothing quality via residual analysis
%   5. Compares nozzle length performance metrics

clear; close all; clc;

% PARAMETERS
target_time = 5;        % Time point in data for velocity/acceleration analysis
sampling_freq = 30;      % Hz/video frame rate
pixels_per_mm = 2.9;     % Calibration from DLTdv8
bl_mm = 260;             % Boat length in mm
sgol_order = 3;          % Savitzky-Golay polynomial order
sgol_framelen = 31;      % Savitzky-Golay frame length (must be odd)
boat_mass_kg = 1.2;      % Boat mass in kg
g = 9.81;                % Gravity constant (m/s²)
voltage = 9.0;           % Single voltage for this dataset

% Folder info
traj_folder = 'Trajectories';
power_folder = 'PowerReadingsCropped';
raw_power_folder = 'RawPower';
nozzle_lengths = [5, 10, 15, 20, 25, 30];  % mm
num_trials = 5;

% Build conditions list
conditions = [arrayfun(@(x) sprintf('%02dmm', x), nozzle_lengths, 'UniformOutput', false), {'Rigid'}];
condition_labels = [arrayfun(@(x) sprintf('%d mm', x), nozzle_lengths, 'UniformOutput', false), {'Rigid'}];

% Color scheme - generate distinct colors for each condition
colors_map = [
    0.8 0.2 0.2;  % 5mm - Red
    0.9 0.5 0.1;  % 10mm - Orange
    0.9 0.8 0.1;  % 15mm - Yellow
    0.2 0.7 0.3;  % 20mm - Green
    0.2 0.5 0.8;  % 25mm - Blue
    0.5 0.2 0.7;  % 30mm - Purple
    0.3 0.3 0.3;  % Rigid - Gray
];
colors.gray = [0.7 0.7 0.7];
colors.green = [0.0 0.6 0.0];
colors.red = [0.8 0.0 0.0];

%% Load and Process Data
fprintf('Processing trajectory data with Savitzky-Golay smoothing...\n');
fprintf('  Filter order: %d, Frame length: %d\n', sgol_order, sgol_framelen);
fit_params = struct();

% Get trajectory files
files = dir(fullfile(traj_folder, '*xypts.csv'));
files = files(~startsWith({files.name}, '._'));

%% Figure 1: Position Data with Smoothed Curves
num_conditions = length(conditions);
subplot_rows = 3;
subplot_cols = 3;
fig1 = figure('Position', [100, 100, 1800, 1200], 'Color', 'w');
set(fig1, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Arial');

plot_idx = 1;
for c_idx = 1:num_conditions
    cond = conditions{c_idx};

    if plot_idx <= subplot_rows * subplot_cols
        subplot(subplot_rows, subplot_cols, plot_idx);
        fit_params = process_condition(cond, files, ...
            traj_folder, sampling_freq, pixels_per_mm, bl_mm, target_time, ...
            sgol_order, sgol_framelen, colors_map(c_idx,:), colors, fit_params, num_trials);
        plot_idx = plot_idx + 1;
    end
end

sgtitle('ChatSweep Analysis - Position with Savitzky-Golay Smoothing (9V)', 'FontSize', 14, 'FontWeight', 'bold');

%% Extract Performance Metrics
[vel_all, accel_all, vel_std, accel_std, vel_traces, accel_traces] = extract_metrics(fit_params, conditions);

%% Load Power Data and Calculate Cost of Transport
[power_all, power_std, power_traces] = load_power_data(power_folder, conditions, num_trials);

% Calculate Cost of Transport: CoT = P / (m * g * v)
% Convert velocity from BL/s to m/s
vel_all_ms = vel_all * (bl_mm / 1000);  % BL/s → m/s
vel_std_ms = vel_std * (bl_mm / 1000);

% Power from mW to W
power_all_W = power_all / 1000;
power_std_W = power_std / 1000;

% Calculate CoT
cot_all = power_all_W ./ (boat_mass_kg * g * vel_all_ms);

% Propagate error for CoT
cot_std = cot_all .* sqrt((power_std_W ./ power_all_W).^2 + (vel_std_ms ./ vel_all_ms).^2);

%% Figure 2: Power Traces
fig2 = figure('Position', [150, 150, 1800, 1000], 'Color', 'w');
set(fig2, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Arial');

plot_idx = 1;
for c_idx = 1:num_conditions
    if plot_idx <= subplot_rows * subplot_cols
        subplot(subplot_rows, subplot_cols, plot_idx);
        hold on;

        % Plot all trials for this condition with error band
        if isfield(power_traces, make_valid_fieldname(conditions{c_idx}))
            trials = power_traces.(make_valid_fieldname(conditions{c_idx}));

            if ~isempty(trials)
                % Find the minimum length across all trials
                min_len = min(cellfun(@length, trials));

                % Truncate all trials to the same length and create matrix
                power_matrix = nan(length(trials), min_len);
                for t = 1:length(trials)
                    if ~isempty(trials{t})
                        power_matrix(t, :) = trials{t}(1:min_len);
                    end
                end

                % Calculate mean and std across trials
                power_mean = mean(power_matrix, 1, 'omitnan');
                power_std_trace = std(power_matrix, 0, 1, 'omitnan');
                samples = 1:min_len;

                % Plot mean line
                plot(samples, power_mean, '-', 'LineWidth', 2, 'Color', colors_map(c_idx,:));

                % Add error band
                fill([samples, fliplr(samples)], ...
                     [power_mean + power_std_trace, fliplr(power_mean - power_std_trace)], ...
                     colors_map(c_idx,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            end

            xlabel('Sample Number', 'FontWeight', 'bold');
            ylabel('Power (mW)', 'FontWeight', 'bold');
            title(condition_labels{c_idx}, 'FontSize', 11);
            box on;
            ylim_vals = ylim;
            ylim([max(0, ylim_vals(1)), ylim_vals(2)]);
        end

        plot_idx = plot_idx + 1;
    end
end

sgtitle('Power Traces - All Trials (9V)', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 3: Velocity Over Time
fig3 = figure('Position', [200, 200, 1800, 1000], 'Color', 'w');
set(fig3, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Arial');

plot_idx = 1;
for c_idx = 1:num_conditions
    if plot_idx <= subplot_rows * subplot_cols
        subplot(subplot_rows, subplot_cols, plot_idx);
        hold on;

        % Plot velocity traces for this condition
        if isfield(vel_traces, make_valid_fieldname(conditions{c_idx}))
            traces = vel_traces.(make_valid_fieldname(conditions{c_idx}));
            if ~isempty(traces) && ~isempty(traces{1})
                time_vec = traces{1}.time;
                vel_vec = traces{1}.velocity;

                % Plot main line
                plot(time_vec, vel_vec, '-', 'LineWidth', 2, 'Color', colors_map(c_idx,:));

                % Add shaded error band using std from target time
                vel_error = vel_std(c_idx);
                if ~isnan(vel_error) && vel_error > 0
                    fill([time_vec; flipud(time_vec)], ...
                         [vel_vec + vel_error; flipud(vel_vec - vel_error)], ...
                         colors_map(c_idx,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end
            end

            xlabel('Time (s)', 'FontWeight', 'bold');
            ylabel('Velocity (BL/s)', 'FontWeight', 'bold');
            title(condition_labels{c_idx}, 'FontSize', 11);
            box on;
            if target_time > 0
                xline(target_time, '--k', 'LineWidth', 1, 'Alpha', 0.5);
            end
        end

        plot_idx = plot_idx + 1;
    end
end

sgtitle('Velocity Over Time - All Trials (9V, S-G Smoothing)', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 4: Acceleration Over Time
fig4 = figure('Position', [250, 250, 1800, 1000], 'Color', 'w');
set(fig4, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Arial');

plot_idx = 1;
for c_idx = 1:num_conditions
    if plot_idx <= subplot_rows * subplot_cols
        subplot(subplot_rows, subplot_cols, plot_idx);
        hold on;

        % Plot acceleration traces for this condition
        if isfield(accel_traces, make_valid_fieldname(conditions{c_idx}))
            traces = accel_traces.(make_valid_fieldname(conditions{c_idx}));
            if ~isempty(traces) && ~isempty(traces{1})
                time_vec = traces{1}.time;
                accel_vec = traces{1}.acceleration;

                % Plot main line
                plot(time_vec, accel_vec, '-', 'LineWidth', 2, 'Color', colors_map(c_idx,:));

                % Add shaded error band using std from target time
                accel_error = accel_std(c_idx);
                if ~isnan(accel_error) && accel_error > 0
                    fill([time_vec; flipud(time_vec)], ...
                         [accel_vec + accel_error; flipud(accel_vec - accel_error)], ...
                         colors_map(c_idx,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end
            end

            yline(0, '--k', 'LineWidth', 1, 'Alpha', 0.5);
            xlabel('Time (s)', 'FontWeight', 'bold');
            ylabel('Acceleration (BL/s²)', 'FontWeight', 'bold');
            title(condition_labels{c_idx}, 'FontSize', 11);
            box on;
            if target_time > 0
                xline(target_time, '--k', 'LineWidth', 1, 'Alpha', 0.5);
            end
        end

        plot_idx = plot_idx + 1;
    end
end

sgtitle('Acceleration Over Time - All Trials (9V, S-G Smoothing)', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 5: Summary Line Plot
fig5 = figure('Position', [300, 300, 1600, 900], 'Color', 'w');
set(fig5, 'DefaultAxesFontSize', 11, 'DefaultAxesFontName', 'Arial');

x_positions = 1:num_conditions;

% Velocity subplot
subplot(2, 2, 1);
hold on;
for i = 1:num_conditions
    if ~isnan(vel_all(i))
        errorbar(x_positions(i), vel_all(i), vel_std(i), 'o', ...
            'MarkerSize', 10, 'MarkerFaceColor', colors_map(i,:), ...
            'MarkerEdgeColor', colors_map(i,:), 'Color', colors_map(i,:), ...
            'LineWidth', 2, 'CapSize', 8);
    end
end
xlabel('Condition', 'FontWeight', 'bold');
ylabel(sprintf('Velocity at t=%.1fs (BL/s)', target_time), 'FontWeight', 'bold');
title('Velocity vs Nozzle Length', 'FontSize', 12);
set(gca, 'XTick', x_positions, 'XTickLabel', condition_labels);
xtickangle(45);
box on;

% Acceleration subplot
subplot(2, 2, 2);
hold on;
for i = 1:num_conditions
    if ~isnan(accel_all(i))
        errorbar(x_positions(i), accel_all(i), accel_std(i), 's', ...
            'MarkerSize', 10, 'MarkerFaceColor', colors_map(i,:), ...
            'MarkerEdgeColor', colors_map(i,:), 'Color', colors_map(i,:), ...
            'LineWidth', 2, 'CapSize', 8);
    end
end
yline(0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Condition', 'FontWeight', 'bold');
ylabel('Acceleration (BL/s^2)', 'FontWeight', 'bold');
title('Acceleration vs Nozzle Length', 'FontSize', 12);
set(gca, 'XTick', x_positions, 'XTickLabel', condition_labels);
xtickangle(45);
box on;

% Power subplot
subplot(2, 2, 3);
hold on;
for i = 1:num_conditions
    if ~isnan(power_all(i))
        errorbar(x_positions(i), power_all(i), power_std(i), '^', ...
            'MarkerSize', 10, 'MarkerFaceColor', colors_map(i,:), ...
            'MarkerEdgeColor', colors_map(i,:), 'Color', colors_map(i,:), ...
            'LineWidth', 2, 'CapSize', 8);
    end
end
xlabel('Condition', 'FontWeight', 'bold');
ylabel('Power (mW)', 'FontWeight', 'bold');
title('Power Consumption vs Nozzle Length', 'FontSize', 12);
set(gca, 'XTick', x_positions, 'XTickLabel', condition_labels);
xtickangle(45);
box on;

% Cost of Transport subplot
subplot(2, 2, 4);
hold on;
for i = 1:num_conditions
    if ~isnan(cot_all(i))
        errorbar(x_positions(i), cot_all(i), cot_std(i), 'd', ...
            'MarkerSize', 10, 'MarkerFaceColor', colors_map(i,:), ...
            'MarkerEdgeColor', colors_map(i,:), 'Color', colors_map(i,:), ...
            'LineWidth', 2, 'CapSize', 8);
    end
end
xlabel('Condition', 'FontWeight', 'bold');
ylabel('Cost of Transport (J/(N·m))', 'FontWeight', 'bold');
title('Cost of Transport vs Nozzle Length', 'FontSize', 12);
set(gca, 'XTick', x_positions, 'XTickLabel', condition_labels);
xtickangle(45);
box on;

sgtitle('Summary: Performance vs Nozzle Length (9V, S-G Smoothing)', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 6: Comparison to Rigid (Ratios)
fig6 = figure('Position', [350, 350, 1600, 900], 'Color', 'w');
set(fig6, 'DefaultAxesFontSize', 11, 'DefaultAxesFontName', 'Arial');

% Get rigid values (last condition)
rigid_idx = num_conditions;
vel_rigid = vel_all(rigid_idx);
accel_rigid = accel_all(rigid_idx);
power_rigid = power_all(rigid_idx);
cot_rigid = cot_all(rigid_idx);

% Calculate ratios (flexible / rigid) with error propagation
vel_ratios = vel_all(1:end-1) / vel_rigid;
accel_ratios = accel_all(1:end-1) / accel_rigid;
power_ratios = power_all(1:end-1) / power_rigid;
cot_ratios = cot_all(1:end-1) / cot_rigid;

% Error propagation for ratios: σ(a/b) = (a/b) * sqrt((σ_a/a)² + (σ_b/b)²)
vel_ratios_std = vel_ratios .* sqrt((vel_std(1:end-1)./vel_all(1:end-1)).^2 + (vel_std(end)/vel_rigid)^2);
accel_ratios_std = accel_ratios .* sqrt((accel_std(1:end-1)./accel_all(1:end-1)).^2 + (accel_std(end)/accel_rigid)^2);
power_ratios_std = power_ratios .* sqrt((power_std(1:end-1)./power_all(1:end-1)).^2 + (power_std(end)/power_rigid)^2);
cot_ratios_std = cot_ratios .* sqrt((cot_std(1:end-1)./cot_all(1:end-1)).^2 + (cot_std(end)/cot_rigid)^2);

x_pos_flex = 1:length(nozzle_lengths);

% Velocity ratio
subplot(2, 2, 1);
plot_ratio_bar(x_pos_flex, vel_ratios, vel_ratios_std, colors_map(1:end-1,:), ...
    nozzle_lengths, 'Velocity Ratio (Flex/Rigid)', colors);

% Acceleration ratio
subplot(2, 2, 2);
plot_ratio_bar(x_pos_flex, accel_ratios, accel_ratios_std, colors_map(1:end-1,:), ...
    nozzle_lengths, 'Acceleration Ratio (Flex/Rigid)', colors);

% Power ratio
subplot(2, 2, 3);
plot_ratio_bar(x_pos_flex, power_ratios, power_ratios_std, colors_map(1:end-1,:), ...
    nozzle_lengths, 'Power Ratio (Flex/Rigid)', colors);

% CoT ratio
subplot(2, 2, 4);
plot_ratio_bar(x_pos_flex, cot_ratios, cot_ratios_std, colors_map(1:end-1,:), ...
    nozzle_lengths, 'CoT Ratio (Flex/Rigid)', colors);

sgtitle('Performance Ratios: Flexible / Rigid (9V, S-G Smoothing)', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 7: Flexible vs Rigid Comparison with Power Insets (Trial 1 only)
fig7 = figure('Position', [400, 400, 1600, 1000], 'Color', 'w');
set(fig7, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Arial');

% Load rigid position data (trial 1)
rigid_file = sprintf('Rigid9VT1Trajxypts.csv');
rigid_data = readmatrix(fullfile(traj_folder, rigid_file));
if isnan(rigid_data(1,1))
    rigid_X = rigid_data(2:end, 1);  % Skip header
else
    rigid_X = rigid_data(:, 1);
end
rigid_X = rigid_X(~isnan(rigid_X));
rigid_X = -rigid_X - min(-rigid_X);
rigid_threshold = find(rigid_X >= 30, 1);
if ~isempty(rigid_threshold)
    rigid_X = movmedian(rigid_X(rigid_threshold:end), 15) / pixels_per_mm;
    rigid_time = (0:length(rigid_X)-1)' / sampling_freq;
end

% Load rigid power data (trial 1)
rigid_power_file = 'Cropped3000_Rigid_T1.txt';
rigid_power = readmatrix(fullfile(power_folder, rigid_power_file));

% Calculate power sampling rate from raw data (trial 1)
rigid_power_sampling_rate = calculate_power_sampling_rate(raw_power_folder, 'Rigid', 1);
if ~isempty(rigid_power) && ~isnan(rigid_power_sampling_rate)
    rigid_power_time = (0:length(rigid_power)-1) / rigid_power_sampling_rate;
else
    rigid_power_time = [];
end

% Plot each flexible nozzle vs rigid
for i = 1:length(nozzle_lengths)
    subplot(3, 2, i);
    hold on;

    % Load flexible position data (trial 1)
    flex_file = sprintf('%02dmm9VT1Trajxypts.csv', nozzle_lengths(i));
    flex_data = readmatrix(fullfile(traj_folder, flex_file));
    if isnan(flex_data(1,1))
        flex_X = flex_data(2:end, 1);  % Skip header
    else
        flex_X = flex_data(:, 1);
    end
    flex_X = flex_X(~isnan(flex_X));
    flex_X = -flex_X - min(-flex_X);
    flex_threshold = find(flex_X >= 30, 1);
    if ~isempty(flex_threshold)
        flex_X = movmedian(flex_X(flex_threshold:end), 15) / pixels_per_mm;
        flex_time = (0:length(flex_X)-1)' / sampling_freq;

        % Plot position traces
        plot(flex_time, flex_X, '-', 'LineWidth', 2, 'Color', colors_map(i,:), ...
            'DisplayName', sprintf('%d mm', nozzle_lengths(i)));
    end

    if exist('rigid_time', 'var')
        plot(rigid_time, rigid_X, '-', 'LineWidth', 2, 'Color', colors_map(end,:), ...
            'DisplayName', 'Rigid');
    end

    xlabel('Time (s)', 'FontWeight', 'bold');
    ylabel('Position (mm)', 'FontWeight', 'bold');
    title(sprintf('%d mm vs Rigid', nozzle_lengths(i)), 'FontSize', 11);
    legend('Location', 'southeast', 'FontSize', 8);
    box on;

    % Create inset axes for power in upper left
    ax_main = gca;
    ax_inset = axes('Position', [ax_main.Position(1) + 0.03, ...
                                  ax_main.Position(2) + ax_main.Position(4) - 0.09, ...
                                  0.10, 0.06]);
    hold(ax_inset, 'on');
    xlim(ax_main,[0,10]);

    % Load flexible power data (trial 1)
    flex_power_file = sprintf('Cropped3000_%dmm_T1.txt', nozzle_lengths(i));
    flex_power = readmatrix(fullfile(power_folder, flex_power_file));

    % Calculate power sampling rate for flexible nozzle
    flex_power_sampling_rate = calculate_power_sampling_rate(raw_power_folder, ...
        sprintf('%02dmm', nozzle_lengths(i)), 1);
    if ~isempty(flex_power) && ~isnan(flex_power_sampling_rate)
        flex_power_time = (0:length(flex_power)-1) / flex_power_sampling_rate;
    else
        flex_power_time = [];
    end

    % Plot power traces in inset
    if ~isempty(flex_power) && ~isempty(flex_power_time)
        plot(ax_inset, flex_power_time, flex_power, '-', 'LineWidth', 1, ...
            'Color', colors_map(i,:));
    end
    if ~isempty(rigid_power) && ~isempty(rigid_power_time)
        plot(ax_inset, rigid_power_time, rigid_power, '-', 'LineWidth', 1, ...
            'Color', colors_map(end,:));
    end

    % Configure inset axes
    box(ax_inset, 'on');
    xlim(ax_inset, [0, 10]);
    ylim(ax_inset, [0, 12000]);
    xticks(ax_inset, 0:2:10);
    yticks(ax_inset, 0:3000:12000);
    set(ax_inset, 'FontSize', 7);
    xlabel(ax_inset, 'Time (s)', 'FontSize', 7);
    ylabel(ax_inset, 'Power (mW)', 'FontSize', 7);
end

sgtitle('Position & Power Comparison: Flexible vs Rigid (Trial 1)', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 8: Performance Ratios vs Inverse C-hat (Point Plot)
fig8 = figure('Position', [450, 450, 1600, 900], 'Color', 'w');
set(fig8, 'DefaultAxesFontSize', 11, 'DefaultAxesFontName', 'Arial');

% Calculate inverse C-hat values for each nozzle length
c_hat_values = zeros(1, length(nozzle_lengths));
for i = 1:length(nozzle_lengths)
    height = nozzle_lengths(i);  % mm
    c_hat = sqrt(270000*0.7/1000/(1000*7/1000))/(height/1000)*0.0125;
    c_hat_values(i) = 1/c_hat;  % Inverse C-hat
end

% Velocity ratio subplot
subplot(2, 2, 1);
plot_ratio_points(c_hat_values, vel_ratios, vel_ratios_std, colors_map(1:end-1,:), ...
    'Velocity Ratio (Flex/Rigid)', colors);

% Acceleration ratio subplot
subplot(2, 2, 2);
plot_ratio_points(c_hat_values, accel_ratios, accel_ratios_std, colors_map(1:end-1,:), ...
    'Acceleration Ratio (Flex/Rigid)', colors);

% Power ratio subplot
subplot(2, 2, 3);
plot_ratio_points(c_hat_values, power_ratios, power_ratios_std, colors_map(1:end-1,:), ...
    'Power Ratio (Flex/Rigid)', colors);

% CoT ratio subplot
subplot(2, 2, 4);
plot_ratio_points(c_hat_values, cot_ratios, cot_ratios_std, colors_map(1:end-1,:), ...
    'CoT Ratio (Flex/Rigid)', colors);

sgtitle('Performance Ratios vs Inverse C-hat: Flexible / Rigid (9V, S-G Smoothing)', 'FontSize', 14, 'FontWeight', 'bold');

%% Print Results
print_results(conditions, condition_labels, vel_all, accel_all, ...
    power_all, cot_all, vel_std, accel_std, power_std, cot_std, target_time);

% Print smoothing quality summary
fprintf('\n========================================\n');
fprintf('SMOOTHING QUALITY (NRMSE values)\n');
fprintf('========================================\n');
for c_idx = 1:num_conditions
    field = make_valid_fieldname(conditions{c_idx});

    if isfield(fit_params, field) && isfield(fit_params.(field), 'nrmse')
        nrmse = fit_params.(field).nrmse;
        fprintf('%s: NRMSE = %.4f', condition_labels{c_idx}, nrmse);
        if nrmse > 0.05
            fprintf(' ⚠');
        end
        fprintf('\n');
    end
end
fprintf('\nNote: NRMSE < 0.05 indicates good smoothing quality\n');

fprintf('\n✓ Analysis complete! Generated 8 figures.\n');

%% ==================== HELPER FUNCTIONS ====================

function fit_params = process_condition(condition, files, ...
    traj_folder, sampling_freq, pixels_per_mm, bl_mm, target_time, ...
    sgol_order, sgol_framelen, color, colors, fit_params, num_trials)
    % Process trajectory data for one condition using Savitzky-Golay smoothing

    % Format condition for filename matching
    if strcmp(condition, 'Rigid')
        prefix = 'Rigid9V';
    else
        prefix = sprintf('%s9V', condition);
    end

    trial_velocities = [];  % Store velocity at target_time for each trial
    trial_accels = [];      % Store acceleration at target_time for each trial
    smoothed_trials = {};   % Store smoothed position data from each trial
    trial_times = {};       % Store time vectors from each trial
    trial_nrmse = [];       % Store NRMSE for each trial

    dt = 1 / sampling_freq;  % Time step

    % Load all trials
    for trial = 1:num_trials
        file_pattern = sprintf('%sT%d', prefix, trial);
        file_idx = find(contains({files.name}, file_pattern), 1);

        if ~isempty(file_idx)
            data = readmatrix(fullfile(traj_folder, files(file_idx).name));

            % Check if data has header
            if isnan(data(1,1))
                X = data(2:end, 1);  % Skip header
            else
                X = data(:, 1);
            end

            X = X(~isnan(X));

            if ~isempty(X)
                % Process trajectory
                X = -X - min(-X);
                threshold_idx = find(X >= 30, 1);
                if ~isempty(threshold_idx)
                    X_raw = movmedian(X(threshold_idx:end), 15) / pixels_per_mm;
                    time = (0:length(X_raw)-1)' / sampling_freq;

                    % Plot raw trial data
                    plot(time, X_raw, '-', 'LineWidth', 1, 'Color', colors.gray);
                    hold on;

                    % Apply Savitzky-Golay smoothing to individual trial
                    if length(X_raw) >= sgol_framelen
                        X_smooth = sgolayfilt(X_raw, sgol_order, sgol_framelen);

                        % Calculate NRMSE for this trial
                        rmse = sqrt(mean((X_raw - X_smooth).^2));
                        nrmse = rmse / (max(X_raw) - min(X_raw));
                        trial_nrmse = [trial_nrmse; nrmse];

                        % Store smoothed data for averaging
                        smoothed_trials{end+1} = X_smooth;
                        trial_times{end+1} = time;

                        % Calculate velocity using gradient (numerical differentiation)
                        velocity = gradient(X_smooth, dt);  % mm/s

                        % Calculate acceleration using gradient of velocity
                        acceleration = gradient(velocity, dt);  % mm/s²

                        % Extract values at target_time
                        if max(time) >= target_time
                            % Find index closest to target_time
                            [~, target_idx] = min(abs(time - target_time));

                            % Store velocity and acceleration at target time
                            v_trial = velocity(target_idx) / bl_mm;  % BL/s
                            a_trial = acceleration(target_idx) / bl_mm;  % BL/s²

                            trial_velocities = [trial_velocities; v_trial];
                            trial_accels = [trial_accels; a_trial];
                        end
                    end
                end
            end
        end
    end

    % Average smoothed curves for visualization
    if ~isempty(smoothed_trials)
        % Find minimum length across all trials
        min_len = min(cellfun(@length, smoothed_trials));

        % Truncate all trials to same length and compute average
        smoothed_matrix = zeros(length(smoothed_trials), min_len);
        for i = 1:length(smoothed_trials)
            smoothed_matrix(i, :) = smoothed_trials{i}(1:min_len);
        end
        avg_smooth = mean(smoothed_matrix, 1);
        avg_time = (0:min_len-1)' / sampling_freq;

        % Average NRMSE across trials
        avg_nrmse = mean(trial_nrmse);

        % Plot averaged smoothed curve
        plot(avg_time, avg_smooth, '-', 'LineWidth', 2.5, 'Color', color);

        % Calculate velocity and acceleration traces from averaged smoothed data
        velocity_smooth = gradient(avg_smooth, dt) / bl_mm;  % BL/s
        acceleration_smooth = gradient(velocity_smooth * bl_mm, dt) / bl_mm;  % BL/s²

        % Store results
        field_name = make_valid_fieldname(condition);
        fit_params.(field_name).nrmse = avg_nrmse;
        fit_params.(field_name).time_fit = avg_time';
        fit_params.(field_name).velocity_fit = velocity_smooth;
        fit_params.(field_name).acceleration_fit = acceleration_smooth;

        % Store metrics from individual trials (averaged)
        if length(trial_velocities) >= 1
            v_at_target_bl = mean(trial_velocities);
            fit_params.(field_name).v_at_target_bl = v_at_target_bl;

            if length(trial_accels) >= 1
                accel_bl = mean(trial_accels);
                fit_params.(field_name).accel_bl = accel_bl;
            else
                fit_params.(field_name).accel_bl = 0;
                accel_bl = 0;
            end

            % Standard deviations from individual trials
            if length(trial_velocities) >= 2
                fit_params.(field_name).v_at_target_std_bl = std(trial_velocities);
            else
                fit_params.(field_name).v_at_target_std_bl = 0;
            end

            if length(trial_accels) >= 2
                fit_params.(field_name).accel_std_bl = std(trial_accels);
            else
                fit_params.(field_name).accel_std_bl = 0;
            end

            % Add annotation with smoothing quality and std dev
            if length(trial_velocities) >= 2
                text(0.98, 0.98, sprintf('v: %.2f±%.2f BL/s\na: %.3f±%.3f BL/s²\nNRMSE: %.4f', ...
                    v_at_target_bl, fit_params.(field_name).v_at_target_std_bl, ...
                    accel_bl, fit_params.(field_name).accel_std_bl, avg_nrmse), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', 'BackgroundColor', 'w', ...
                    'EdgeColor', 'k', 'FontSize', 9, 'Margin', 3);
            else
                text(0.98, 0.98, sprintf('v: %.2f BL/s\na: %.3f BL/s²\nNRMSE: %.4f', ...
                    v_at_target_bl, accel_bl, avg_nrmse), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', 'BackgroundColor', 'w', ...
                    'EdgeColor', 'k', 'FontSize', 9, 'Margin', 3);
            end
        end
    end

    % Format subplot
    xlabel('Time (s)', 'FontWeight', 'bold');
    ylabel('Position (mm)', 'FontWeight', 'bold');
    title(strrep(condition, '_', ' '), 'FontSize', 11);
    box on;
    if ~isempty(smoothed_trials)
        xlim([0 max(avg_time)]);
    end
end

function [vel_all, accel_all, vel_std, accel_std, vel_traces, accel_traces] = extract_metrics(fit_params, conditions)
    % Extract velocity and acceleration arrays from fit parameters
    % Also returns time series traces

    num_cond = length(conditions);
    vel_all = nan(1, num_cond);
    accel_all = nan(1, num_cond);
    vel_std = nan(1, num_cond);
    accel_std = nan(1, num_cond);
    vel_traces = struct();
    accel_traces = struct();

    for i = 1:num_cond
        field = make_valid_fieldname(conditions{i});

        if isfield(fit_params, field)
            if isfield(fit_params.(field), 'v_at_target_bl')
                vel_all(i) = fit_params.(field).v_at_target_bl;
            end
            if isfield(fit_params.(field), 'accel_bl')
                accel_all(i) = fit_params.(field).accel_bl;
            end
            if isfield(fit_params.(field), 'v_at_target_std_bl')
                vel_std(i) = fit_params.(field).v_at_target_std_bl;
            end
            if isfield(fit_params.(field), 'accel_std_bl')
                accel_std(i) = fit_params.(field).accel_std_bl;
            end

            % Extract time series traces
            if isfield(fit_params.(field), 'time_fit') && isfield(fit_params.(field), 'velocity_fit')
                vel_traces.(field) = {struct('time', fit_params.(field).time_fit, ...
                                              'velocity', fit_params.(field).velocity_fit)};
            end
            if isfield(fit_params.(field), 'time_fit') && isfield(fit_params.(field), 'acceleration_fit')
                accel_traces.(field) = {struct('time', fit_params.(field).time_fit, ...
                                                'acceleration', fit_params.(field).acceleration_fit)};
            end
        end
    end
end

function [power_all, power_std, power_traces] = load_power_data(power_folder, conditions, num_trials)
    % Load and average power consumption data for each condition
    % Also returns raw power traces for plotting

    num_cond = length(conditions);
    power_all = nan(1, num_cond);
    power_std = nan(1, num_cond);
    power_traces = struct();

    for c_idx = 1:num_cond
        condition = conditions{c_idx};

        % Load all trials
        trial_powers = [];
        trial_traces = {};

        for trial = 1:num_trials
            % Construct filename for cropped data
            if strcmp(condition, 'Rigid')
                filename = sprintf('Cropped3000_Rigid_T%d.txt', trial);
            else
                % Extract number from condition (e.g., '05mm' -> '5')
                nozzle_num = str2double(condition(1:2));
                filename = sprintf('Cropped3000_%dmm_T%d.txt', nozzle_num, trial);
            end
            filepath = fullfile(power_folder, filename);

            if isfile(filepath)
                try
                    % Read cropped power file (single row, comma-separated)
                    data = readmatrix(filepath);

                    if ~isempty(data)
                        % Data is a row vector of power values in mW
                        % Average power over entire cropped period
                        trial_powers(end+1) = mean(data(:));
                        % Store the trace
                        trial_traces{end+1} = data(:)';
                    end
                catch ME
                    fprintf('Warning: Error reading %s: %s\n', filename, ME.message);
                end
            end
        end

        % Calculate mean and std across trials
        if ~isempty(trial_powers)
            power_all(c_idx) = mean(trial_powers);
            if length(trial_powers) >= 2
                power_std(c_idx) = std(trial_powers);
            else
                power_std(c_idx) = 0;
            end
        end

        % Store traces
        field_name = make_valid_fieldname(condition);
        power_traces.(field_name) = trial_traces;
    end
end

function plot_ratio_bar(x_pos, ratios, ratios_std, colors_arr, nozzle_lengths, ylabel_str, colors)
    % Create bar chart for ratios with color coding and error bars

    b = bar(x_pos, ratios);
    b.FaceColor = 'flat';
    for i = 1:length(x_pos)
        b.CData(i,:) = colors_arr(i,:);
    end
    b.EdgeColor = 'k';
    b.LineWidth = 1.2;

    hold on;

    % Add error bars
    errorbar(x_pos, ratios, ratios_std, 'k.', 'LineWidth', 1.5, 'CapSize', 10, ...
        'HandleVisibility', 'off');

    yline(1, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Add percentage labels (positioned higher)
    for i = 1:length(ratios)
        if ~isnan(ratios(i))
            pct = (ratios(i) - 1) * 100;
            % Increase y_offset to position labels higher
            y_range = max(ratios) - min(ratios);
            y_offset = max(0.08, y_range * 0.08);
            if ratios(i) < 1
                label_color = colors.red;
            else
                label_color = colors.green;
            end
            text(x_pos(i), ratios(i) + y_offset, sprintf('%+.1f%%', pct), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
                'FontSize', 10, 'Color', label_color);
        end
    end

    xlabel('Nozzle Length (mm)', 'FontWeight', 'bold');
    ylabel(ylabel_str, 'FontWeight', 'bold');
    set(gca, 'XTick', x_pos, 'XTickLabel', arrayfun(@num2str, nozzle_lengths, 'UniformOutput', false));
    box on;

    % Adjust ylim to accommodate labels
    ylims = ylim;
    ylim([ylims(1), ylims(2) * 1.12]);
end

function plot_ratio_points(c_hat_values, ratios, ratios_std, colors_arr, ylabel_str, colors)
    % Create scatter plot for ratios with blue points, lines, and error bars

    hold on;

    % Plot horizontal reference line at y=1
    yline(1, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Define blue color
    blue_color = [0.2 0.5 0.8];

    % Plot connecting line
    plot(c_hat_values, ratios, '-', 'Color', blue_color, 'LineWidth', 2);

    % Plot each point with error bar
    for i = 1:length(c_hat_values)
        % Plot error bar
        errorbar(c_hat_values(i), ratios(i), ratios_std(i), 'k', ...
            'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');

        % Plot point (all blue)
        scatter(c_hat_values(i), ratios(i), 100, blue_color, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2);

        % Add percentage labels
        if ~isnan(ratios(i))
            pct = (ratios(i) - 1) * 100;
            y_range = max(ratios) - min(ratios);
            y_offset = max(0.08, y_range * 0.08);
            if ratios(i) < 1
                label_color = colors.red;
            else
                label_color = colors.green;
            end
            text(c_hat_values(i), ratios(i) + y_offset, sprintf('%+.1f%%', pct), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
                'FontSize', 10, 'Color', label_color);
        end
    end

    xlabel('Inverse C-hat (1/\^C)', 'FontWeight', 'bold');
    ylabel(ylabel_str, 'FontWeight', 'bold');
    box on;
    grid on;

    % Adjust ylim to accommodate labels
    ylims = ylim;
    ylim([ylims(1), ylims(2) * 1.12]);

    % Adjust xlim for better visualization
    x_range = max(c_hat_values) - min(c_hat_values);
    xlim([min(c_hat_values) - 0.1*x_range, max(c_hat_values) + 0.1*x_range]);

    % Set aspect ratio to [1.5, 1, 1]
    pbaspect([1.5 1 1]);
end

function print_results(conditions, condition_labels, vel_all, accel_all, ...
    power_all, cot_all, vel_std, accel_std, power_std, cot_std, target_time)
    % Print analysis results to console

    fprintf('\n========================================\n');
    fprintf('PERFORMANCE SUMMARY\n');
    fprintf('========================================\n');
    fprintf('Analysis time point: t = %.1fs\n', target_time);
    fprintf('Voltage: 9.0V\n\n');

    fprintf('%-12s | %-15s | %-15s | %-15s | %-15s\n', ...
        'Condition', 'Vel (BL/s)', 'Accel (BL/s²)', 'Power (mW)', 'CoT (J/N·m)');
    fprintf('-------------|-----------------|-----------------|-----------------|----------------\n');

    for i = 1:length(conditions)
        if ~isnan(vel_all(i))
            fprintf('%-12s | %7.2f ± %5.2f | %7.3f ± %5.3f | %7.1f ± %5.1f | %7.3f ± %5.3f\n', ...
                condition_labels{i}, vel_all(i), vel_std(i), accel_all(i), accel_std(i), ...
                power_all(i), power_std(i), cot_all(i), cot_std(i));
        end
    end
    fprintf('\n');

    % Print comparison to rigid
    rigid_idx = length(conditions);
    fprintf('Comparison to Rigid:\n');
    fprintf('%-12s | %-10s | %-10s | %-10s | %-10s\n', ...
        'Condition', 'Vel Δ', 'Accel Δ', 'Power Δ', 'CoT Δ');
    fprintf('-------------|-----------|-----------|-----------|------------\n');

    for i = 1:rigid_idx-1
        if ~isnan(vel_all(i))
            pct_vel = ((vel_all(i) - vel_all(rigid_idx)) / vel_all(rigid_idx)) * 100;
            pct_accel = ((accel_all(i) - accel_all(rigid_idx)) / accel_all(rigid_idx)) * 100;
            pct_power = ((power_all(i) - power_all(rigid_idx)) / power_all(rigid_idx)) * 100;
            pct_cot = ((cot_all(i) - cot_all(rigid_idx)) / cot_all(rigid_idx)) * 100;

            fprintf('%-12s | %+9.1f%% | %+9.1f%% | %+9.1f%% | %+9.1f%%\n', ...
                condition_labels{i}, pct_vel, pct_accel, pct_power, pct_cot);
        end
    end
    fprintf('\n');
end

function field_name = make_valid_fieldname(condition)
    % Convert condition name to valid MATLAB field name
    % Field names must start with a letter, so prepend 'n' to numeric starts

    field_name = strrep(condition, '.', '_');

    % If starts with a digit, prepend 'n'
    if ~isempty(field_name) && isstrprop(field_name(1), 'digit')
        field_name = ['ex' field_name];
    end
end

function sampling_rate = calculate_power_sampling_rate(raw_power_folder, condition, trial)
    % Calculate power sampling rate from raw power data timestamps
    % Returns sampling rate in Hz, or NaN if calculation fails

    % Construct raw power filename
    if strcmp(condition, 'Rigid')
        raw_filename = sprintf('Rigid9VT%d.txt', trial);
    else
        raw_filename = sprintf('%s9VT%d.txt', condition, trial);
    end
    raw_filepath = fullfile(raw_power_folder, raw_filename);

    % Default return value
    sampling_rate = NaN;

    if ~isfile(raw_filepath)
        fprintf('Warning: Raw power file not found: %s\n', raw_filename);
        return;
    end

    try
        % Read the file as text
        fid = fopen(raw_filepath, 'r');
        if fid == -1
            fprintf('Warning: Could not open file: %s\n', raw_filename);
            return;
        end

        % Read first line to get initial timestamp
        first_line = fgetl(fid);
        if ischar(first_line) && ~isempty(first_line)
            % Extract timestamp from first line (format: HH:MM:SS.mmm -> ...)
            first_timestamp_str = regexp(first_line, '(\d+):(\d+):(\d+\.\d+)', 'tokens', 'once');
            if ~isempty(first_timestamp_str)
                first_hours = str2double(first_timestamp_str{1});
                first_mins = str2double(first_timestamp_str{2});
                first_secs = str2double(first_timestamp_str{3});
                first_time = first_hours * 3600 + first_mins * 60 + first_secs;
            else
                fclose(fid);
                fprintf('Warning: Could not parse first timestamp in: %s\n', raw_filename);
                return;
            end
        else
            fclose(fid);
            fprintf('Warning: Empty or invalid first line in: %s\n', raw_filename);
            return;
        end

        % Count total lines and get last line
        line_count = 1;
        last_line = first_line;
        while ~feof(fid)
            current_line = fgetl(fid);
            if ischar(current_line) && ~isempty(current_line)
                last_line = current_line;
                line_count = line_count + 1;
            end
        end
        fclose(fid);

        % Extract timestamp from last line
        last_timestamp_str = regexp(last_line, '(\d+):(\d+):(\d+\.\d+)', 'tokens', 'once');
        if ~isempty(last_timestamp_str)
            last_hours = str2double(last_timestamp_str{1});
            last_mins = str2double(last_timestamp_str{2});
            last_secs = str2double(last_timestamp_str{3});
            last_time = last_hours * 3600 + last_mins * 60 + last_secs;
        else
            fprintf('Warning: Could not parse last timestamp in: %s\n', raw_filename);
            return;
        end

        % Calculate sampling rate
        time_duration = last_time - first_time;
        if time_duration > 0 && line_count > 1
            sampling_rate = (line_count - 1) / time_duration;  % Hz
            fprintf('Power sampling rate for %s T%d: %.2f Hz (%.1f s, %d samples)\n', ...
                condition, trial, sampling_rate, time_duration, line_count);
        else
            fprintf('Warning: Invalid time duration or line count in: %s\n', raw_filename);
        end

    catch ME
        fprintf('Warning: Error calculating power sampling rate for %s: %s\n', ...
            raw_filename, ME.message);
    end
end
