% ChatSweep Data Analysis using Savitzky-Golay smoothing
% Compares different nozzle lengths at 9V

% METHOD:
%   1. Applies Savitzky-Golay smoothing filter to position data
%   2. Computes velocity via gradient of smoothed position
%   3. Computes acceleration via gradient of smoothed velocity
%   4. Reports smoothing parameters (window size, polynomial order)

clear; close all; clc;

% PARAMETERS
target_time = 3;        % Time point in data for velocity/acceleration analysis
video_fps = 240;        % Video recording frame rate (fps)
sampling_freq = 1.2;    % Tracking frequency (Hz) - data tracked once every 200 frames (240/200 = 1.2 Hz)
pixels_per_mm = 2.6;    % Calibration from DLTdv8
bl_mm = 65;             % Boat length in mm
sgol_order = 3;         % Savitzky-Golay polynomial order
sgol_framelen = 11;     % Savitzky-Golay frame length (must be odd, > sgol_order)
boat_mass_kg = 0.1;      % Boat mass in kg
g = 9.81;               % Gravity constant (m/s²)
voltage = 3.7;          % Single voltage for this dataset

% Folder info
traj_folder = '21Nov25TrajectoryData';
power_folder = 'PowerReadingsCropped';  % Not used with current dataset
raw_power_folder = 'RawPower';  % Not used with current dataset
trials = [1, 2, 4];  % Available trials

% Build conditions list - only Flexible Nozzle (FN) and Rigid Nozzle (RN)
conditions = {'FN', 'RN'};
condition_labels = {'Flexible Nozzle', 'Rigid Nozzle'};

% Color scheme - colors for each condition
colors_map = [
    0.2 0.5 0.8;  % FN - Blue
    0.3 0.3 0.3;  % RN - Gray
];
colors.gray = [0.7 0.7 0.7];
colors.green = [0.0 0.6 0.0];
colors.red = [0.8 0.0 0.0];

%% Load and Process Data
fprintf('Processing trajectory data (Savitzky-Golay smoothing)...\n');
fit_params = struct();

% Get trajectory files
files = dir(fullfile(traj_folder, '*xypts.tsv'));
files = files(~startsWith({files.name}, '._'));

%% Figure 1: Position Data with Smoothed Curves
num_conditions = length(conditions);
subplot_rows = 1;
subplot_cols = 2;
fig1 = figure('Position', [100, 100, 1200, 500], 'Color', 'w');
set(fig1, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Arial');

plot_idx = 1;
for c_idx = 1:num_conditions
    cond = conditions{c_idx};

    subplot(subplot_rows, subplot_cols, plot_idx);
    fit_params = process_condition(cond, files, trials, ...
        traj_folder, video_fps, pixels_per_mm, bl_mm, target_time, ...
        sgol_order, sgol_framelen, colors_map(c_idx,:), colors, fit_params);
    plot_idx = plot_idx + 1;
end

sgtitle('ChatSweep Analysis - Position with Savitzky-Golay Smoothing (9V)', 'FontSize', 14, 'FontWeight', 'bold');

%% Extract Performance Metrics
[vel_all, accel_all, vel_std, accel_std, vel_traces, accel_traces] = extract_metrics(fit_params, conditions);

%% Power Data Not Available
% Commenting out power and CoT analysis - no power data in current dataset
% [power_all, power_std, power_traces] = load_power_data(power_folder, conditions, trials);
power_all = nan(1, num_conditions);
power_std = nan(1, num_conditions);
power_traces = struct();
cot_all = nan(1, num_conditions);
cot_std = nan(1, num_conditions);

%% Figure 2: Velocity Over Time
fig2 = figure('Position', [200, 200, 1200, 500], 'Color', 'w');
set(fig2, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Arial');

plot_idx = 1;
for c_idx = 1:num_conditions
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

sgtitle('Velocity Over Time - Savitzky-Golay Smoothing (9V)', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 3: Acceleration Over Time
fig3 = figure('Position', [250, 250, 1200, 500], 'Color', 'w');
set(fig3, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Arial');

plot_idx = 1;
for c_idx = 1:num_conditions
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

sgtitle('Acceleration Over Time - Savitzky-Golay Smoothing (9V)', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 4: Summary Comparison
fig4 = figure('Position', [300, 300, 1200, 500], 'Color', 'w');
set(fig4, 'DefaultAxesFontSize', 11, 'DefaultAxesFontName', 'Arial');

x_positions = 1:num_conditions;

% Velocity subplot
subplot(1, 2, 1);
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
title('Velocity Comparison', 'FontSize', 12);
set(gca, 'XTick', x_positions, 'XTickLabel', condition_labels);
xtickangle(45);
box on;

% Acceleration subplot
subplot(1, 2, 2);
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
ylabel('Acceleration (BL/s²)', 'FontWeight', 'bold');
title('Acceleration Comparison', 'FontSize', 12);
set(gca, 'XTick', x_positions, 'XTickLabel', condition_labels);
xtickangle(45);
box on;

sgtitle('Summary: Flexible vs Rigid Nozzle Performance (9V)', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 5: Velocity Bar Chart with Percent Increase
fig5 = figure('Position', [350, 350, 800, 600], 'Color', 'w');
set(fig5, 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Arial');

% Get velocities
vel_fn = vel_all(1);  % Flexible Nozzle
vel_rn = vel_all(2);  % Rigid Nozzle
vel_fn_std = vel_std(1);
vel_rn_std = vel_std(2);

% Calculate percent increase
pct_increase = ((vel_fn - vel_rn) / vel_rn) * 100;

% Create bar chart
x_pos = [1, 2];
bar_data = [vel_fn, vel_rn];
bar_colors = [colors_map(1,:); colors_map(2,:)];

b = bar(x_pos, bar_data);
b.FaceColor = 'flat';
b.CData = bar_colors;
b.EdgeColor = 'k';
b.LineWidth = 1.5;

hold on;

% Add error bars
errorbar(x_pos, bar_data, [vel_fn_std, vel_rn_std], 'k.', 'LineWidth', 2, 'CapSize', 15);

% Add value labels on top of bars
text(1, vel_fn + vel_fn_std + 0.05*max(bar_data), sprintf('%.2f BL/s', vel_fn), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
text(2, vel_rn + vel_rn_std + 0.05*max(bar_data), sprintf('%.2f BL/s', vel_rn), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);

% Add percent increase annotation
if pct_increase > 0
    annotation_color = colors.green;
    arrow_direction = '\uparrow';
else
    annotation_color = colors.red;
    arrow_direction = '\downarrow';
end
text(1.5, max(bar_data)*0.85, sprintf('%s %+.1f%%', arrow_direction, pct_increase), ...
    'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold', ...
    'Color', annotation_color, 'BackgroundColor', [1 1 1 0.8], 'EdgeColor', annotation_color, ...
    'LineWidth', 2, 'Margin', 5);

% Format plot
ylabel('Velocity at t=3s (BL/s)', 'FontWeight', 'bold', 'FontSize', 13);
set(gca, 'XTick', x_pos, 'XTickLabel', {'Flexible Nozzle', 'Rigid Nozzle'});
xtickangle(0);
ylim([0, 0.54]);%max(bar_data)*1.25]);
box on;
grid on;
grid(gca, 'minor');
title('Velocity Comparison: Flexible vs Rigid Nozzle', 'FontSize', 14, 'FontWeight', 'bold');

%% Figures 6-8: Skipped for 2-condition dataset
% The following figures (ratio comparisons, C-hat analysis) are designed for
% multiple nozzle lengths and have been commented out for the FN vs RN comparison

% Figures 6-8 commented out for FN vs RN dataset

%% Print Results
print_results(conditions, condition_labels, vel_all, accel_all, ...
    power_all, cot_all, vel_std, accel_std, power_std, cot_std, target_time);

% Print smoothing parameters summary
fprintf('\n========================================\n');
fprintf('SMOOTHING PARAMETERS\n');
fprintf('========================================\n');
fprintf('Method: Savitzky-Golay filter\n');
fprintf('Polynomial order: %d\n', sgol_order);
fprintf('Frame length: %d\n', sgol_framelen);
fprintf('\n');

for c_idx = 1:num_conditions
    field = make_valid_fieldname(conditions{c_idx});

    if isfield(fit_params, field) && isfield(fit_params.(field), 'residual_std')
        fprintf('%s:\n', condition_labels{c_idx});

        % Print individual trial residual std values
        if isfield(fit_params.(field), 'trial_residual_std')
            trial_res = fit_params.(field).trial_residual_std;
            for i = 1:length(trial_res)
                fprintf('  Trial %d: residual std = %.3f mm\n', trials(i), trial_res(i));
            end
        end

        % Print average residual std
        fprintf('  Average residual std: %.3f mm\n\n', fit_params.(field).residual_std);
    end
end

fprintf('\n✓ Analysis complete! Generated 5 figures.\n');

%% ==================== HELPER FUNCTIONS ====================

function fit_params = process_condition(condition, files, trials, ...
    traj_folder, video_fps, pixels_per_mm, bl_mm, target_time, ...
    sgol_order, sgol_framelen, color, colors, fit_params)
    % Process trajectory data for one condition using Savitzky-Golay smoothing
    % video_fps: Video recording frame rate (fps) - used to convert frame numbers to time

    % Format condition for filename matching (FN.T# or RN.T#)
    if strcmp(condition, 'RN')
        prefix = 'RN.T';
    else
        prefix = 'FN.T';
    end

    all_time = [];
    all_position = [];
    trial_velocities = [];  % Store velocity at target_time for each trial
    trial_accels = [];      % Store acceleration at target_time for each trial
    trial_residual_std = []; % Store residual std for each trial

    % Load all trials
    for trial = trials
        file_pattern = sprintf('%s%d', prefix, trial);
        file_idx = find(contains({files.name}, file_pattern), 1);

        if ~isempty(file_idx)
            % Read TSV file - sparse format with frame, dimension, value
            data = readmatrix(fullfile(traj_folder, files(file_idx).name), ...
                'FileType', 'text', 'Delimiter', '\t');

            % Parse sparse format: col1=frame, col2=dimension (1=X, 2=Y), col3=value
            % Extract X coordinates (dimension = 1)
            x_mask = data(:, 2) == 1;
            frames = data(x_mask, 1);
            X = data(x_mask, 3);

            X = X(~isnan(X));
            frames = frames(~isnan(frames));

            if ~isempty(X)
                % Convert frame numbers to time using video frame rate
                time = frames / video_fps;

                % Normalize time to start at 0
                time = time - min(time);

                % Process trajectory (convert pixel coordinates)
                X = X / pixels_per_mm;

                % Normalize position to start at 0
                X = X - X(1);

                % Flip X axis so velocity is positive (moving in positive direction)
                X = -X;

                % Sort by time
                [time, sort_idx] = sort(time);
                X = X(sort_idx);

                % Apply Savitzky-Golay smoothing to individual trial
                if length(X) >= sgol_framelen && max(time) >= target_time
                    % Smooth position
                    X_smooth = sgolayfilt(X, sgol_order, sgol_framelen);

                    % Calculate residual std for this trial
                    res_std = std(X - X_smooth);
                    trial_residual_std = [trial_residual_std; res_std];

                    % Calculate dt (assume uniform spacing)
                    dt = mean(diff(time));

                    % Compute velocity via gradient of smoothed position
                    v_trial = gradient(X_smooth, dt);

                    % Compute acceleration via gradient of velocity
                    a_trial = gradient(v_trial, dt);

                    % Interpolate to get values at target_time
                    v_at_target = interp1(time, v_trial, target_time, 'linear');
                    a_at_target = interp1(time, a_trial, target_time, 'linear');

                    % Convert to BL/s
                    trial_velocities = [trial_velocities; v_at_target / bl_mm];
                    trial_accels = [trial_accels; a_at_target / bl_mm];
                end

                % Plot trial data
                plot(time, X, '-', 'LineWidth', 1, 'Color', colors.gray);
                hold on;

                all_time = [all_time; time];
                all_position = [all_position; X];
            end
        end
    end

    % Apply Savitzky-Golay smoothing to combined data
    if ~isempty(all_time)
        [all_time, idx] = sort(all_time);
        all_position = all_position(idx);

        % Ensure frame length does not exceed data length, and is odd
        actual_framelen = sgol_framelen;
        if length(all_position) < actual_framelen
            actual_framelen = length(all_position);
            if mod(actual_framelen, 2) == 0
                actual_framelen = actual_framelen - 1;
            end
        end
        if actual_framelen <= sgol_order
            actual_framelen = sgol_order + 2;
            if mod(actual_framelen, 2) == 0
                actual_framelen = actual_framelen + 1;
            end
        end

        % Apply Savitzky-Golay filter to position
        x_smooth = sgolayfilt(all_position, sgol_order, actual_framelen);

        % Calculate residual std for combined data
        residual_std_combined = std(all_position - x_smooth);

        % Calculate dt (average time step)
        dt = mean(diff(all_time));

        % Compute velocity and acceleration via numerical differentiation of smoothed data
        v_smooth = gradient(x_smooth, dt);  % mm/s
        a_smooth = gradient(v_smooth, dt);  % mm/s²

        % Convert to BL/s
        v_smooth_bl = v_smooth / bl_mm;
        a_smooth_bl = a_smooth / bl_mm;

        % Plot smoothed curve
        plot(all_time, x_smooth, '-', 'LineWidth', 2.5, 'Color', color);

        % Store results
        field_name = make_valid_fieldname(condition);
        fit_params.(field_name).residual_std = residual_std_combined;
        fit_params.(field_name).trial_residual_std = trial_residual_std;
        fit_params.(field_name).time_fit = all_time;
        fit_params.(field_name).velocity_fit = v_smooth_bl;
        fit_params.(field_name).acceleration_fit = a_smooth_bl;

        % Use average of individual trial velocities
        if target_time <= max(all_time) && length(trial_velocities) >= 1
            v_at_target_bl = mean(trial_velocities);  % BL/s
            fit_params.(field_name).v_at_target_bl = v_at_target_bl;

            if length(trial_accels) >= 1
                accel_bl = mean(trial_accels);  % BL/s²
                fit_params.(field_name).accel_bl = accel_bl;
            else
                fit_params.(field_name).accel_bl = 0;
                accel_bl = 0;
            end

            % Standard deviations
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

            % Add annotation with residual info
            res_str = '';
            for i = 1:length(trial_residual_std)
                res_str = [res_str sprintf('res(T%d): %.2f mm\n', trials(i), trial_residual_std(i))];
            end

            if length(trial_velocities) >= 2
                text(0.98, 0.98, sprintf('v: %.2f±%.2f BL/s\na: %.3f±%.3f BL/s²\n%s', ...
                    v_at_target_bl, fit_params.(field_name).v_at_target_std_bl, ...
                    accel_bl, fit_params.(field_name).accel_std_bl, res_str), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', 'BackgroundColor', 'w', ...
                    'EdgeColor', 'k', 'FontSize', 9, 'Margin', 3);
            else
                text(0.98, 0.98, sprintf('v: %.2f BL/s\na: %.3f BL/s²\n%s', ...
                    v_at_target_bl, accel_bl, res_str), ...
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
    xlim([0 10]);
    ylim([0 350]);
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
