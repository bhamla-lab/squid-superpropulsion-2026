% Dat analysis using polynomial fitting

% METHOD:
%   1. Normalizes data (center & scale) for numerical stability
%   2. Fits polynomial of adjustable order (default: cubic)
%   3. Calculates derivatives analytically using polyder()
%   4. Transforms results back to original scale
%   5. Reports R² goodness-of-fit for each condition

clear; close all; clc;

% PARAMETERS 
target_time = 5;        % Time point in data for velocity/acceleration analysis
sampling_freq = 30;     % Hz/video frame rate
pixels_per_mm = 2.9;    % Calibration from DLTdv8
bl_mm = 260;            % Boat length in mm
poly_order = 3;         % Polynomial order for fitting
boat_mass_kg = 1.2;     % Boat mass in kg
g = 9.81;               % Gravity constant (m/s²)

% Folder info
traj_folder = 'Trajectories';
power_folder = 'PowerReadingsCropped';
voltages = [7.8, 8.1, 8.4, 8.7, 9.0, 9.3];
frequencies = [16, 17, 18, 21, 22, 24];  % Hz corresponding to each voltage
voltage_labels = arrayfun(@(v) sprintf('%.1fV', v), voltages, 'UniformOutput', false);
frequency_labels = arrayfun(@(f) sprintf('%dHz', f), frequencies, 'UniformOutput', false);
conditions = {'Flex15mm', 'Rigid'};

% Color scheme
colors.flex = [0.2 0.6 0.8];
colors.rigid = [0.8 0.4 0.2];
colors.gray = [0.7 0.7 0.7];
colors.green = [0.0 0.6 0.0];
colors.red = [0.8 0.0 0.0];

%% Load and Process Data
fprintf('Processing trajectory data...\n');
fit_params = struct();

% Get trajectory files
files = dir(fullfile(traj_folder, '*xypts.csv'));
files = files(~startsWith({files.name}, '._'));

%% Figure 1: Position Data with Fitted Curves
fig1 = figure('Position', [100, 100, 1800, 1200], 'Color', 'w');
set(fig1, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Arial');

plot_idx = 1;
for cond = conditions
    for v_idx = 1:length(voltages)
        subplot(3, 4, plot_idx);
        fit_params = process_condition(cond{1}, voltages(v_idx), files, ...
            traj_folder, sampling_freq, pixels_per_mm, bl_mm, target_time, ...
            poly_order, colors, fit_params);
        plot_idx = plot_idx + 1;
    end
end

sgtitle('NoWire Validation - Position with Fitted Curves', 'FontSize', 14, 'FontWeight', 'bold');

%% Extract Performance Metrics
[vel_flex, vel_rigid, accel_flex, accel_rigid, vel_flex_std, vel_rigid_std, accel_flex_std, accel_rigid_std] = extract_metrics(fit_params, voltages);

%% Load Power Data and Calculate Cost of Transport
[power_flex, power_rigid, power_flex_std, power_rigid_std] = load_power_data(power_folder, voltages, conditions);

% Calculate Cost of Transport: CoT = P / (m * g * v)
% P in Watts, m in kg, g in m/s², v in m/s
% Result in dimensionless units (Joules per Newton-meter)

% Convert velocity from BL/s to m/s
vel_flex_ms = vel_flex * (bl_mm / 1000);  % BL/s → m/s
vel_rigid_ms = vel_rigid * (bl_mm / 1000);

% Power from mW to W
power_flex_W = power_flex / 1000;
power_rigid_W = power_rigid / 1000;

% Calculate CoT
cot_flex = power_flex_W ./ (boat_mass_kg * g * vel_flex_ms);
cot_rigid = power_rigid_W ./ (boat_mass_kg * g * vel_rigid_ms);

% Propagate error for CoT using error propagation formula
% For CoT = P / (m * g * v), where m and g are constants:
% CoT = P / (k * v), where k = m * g
% Using δ(f/g) = f/g * sqrt((δf/f)² + (δg/g)²)
% σ_CoT = CoT * sqrt((σ_P/P)² + (σ_v/v)²)

% Convert velocity std from BL/s to m/s
vel_flex_std_ms = vel_flex_std * (bl_mm / 1000);
vel_rigid_std_ms = vel_rigid_std * (bl_mm / 1000);

% Power std from mW to W
power_flex_std_W = power_flex_std / 1000;
power_rigid_std_W = power_rigid_std / 1000;

% Calculate relative errors and propagate
cot_flex_std = cot_flex .* sqrt((power_flex_std_W ./ power_flex_W).^2 + ...
                                 (vel_flex_std_ms ./ vel_flex_ms).^2);
cot_rigid_std = cot_rigid .* sqrt((power_rigid_std_W ./ power_rigid_W).^2 + ...
                                   (vel_rigid_std_ms ./ vel_rigid_ms).^2);

%% Figure 2: Summary Line Plots (4 panels)
fig2 = figure('Position', [150, 150, 1600, 900], 'Color', 'w');
set(fig2, 'DefaultAxesFontSize', 11, 'DefaultAxesFontName', 'Arial');

% Velocity subplot
subplot(2, 2, 1);
plot_comparison_with_error(voltages, vel_flex, vel_rigid, vel_flex_std, vel_rigid_std, colors, ...
    sprintf('Velocity at t=%.1fs (BL/s)', target_time), ...
    'Velocity vs Voltage');

% Acceleration subplot
subplot(2, 2, 2);
plot_comparison_with_error(voltages, accel_flex, accel_rigid, accel_flex_std, accel_rigid_std, colors, ...
    'Acceleration (BL/s^2)', 'Acceleration vs Voltage');
yline(0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');

% Power subplot
subplot(2, 2, 3);
plot_comparison_with_error(voltages, power_flex, power_rigid, power_flex_std, power_rigid_std, colors, ...
    'Power (mW)', 'Power Consumption vs Voltage');

% Cost of Transport subplot
subplot(2, 2, 4);
plot_comparison_with_error(voltages, cot_flex, cot_rigid, cot_flex_std, cot_rigid_std, colors, ...
    'Cost of Transport (J/(N·m))', 'Cost of Transport vs Voltage');

sgtitle('Summary: Performance vs Voltage', 'FontSize', 14, 'FontWeight', 'bold');

%% Figures 3-8: Individual Voltage Bar Charts
for v_idx = 1:length(voltages)
    if ~isnan(vel_flex(v_idx)) && ~isnan(vel_rigid(v_idx))
        create_voltage_figure(v_idx, voltage_labels{v_idx}, ...
            vel_flex(v_idx), vel_rigid(v_idx), ...
            accel_flex(v_idx), accel_rigid(v_idx), ...
            power_flex(v_idx), power_rigid(v_idx), ...
            cot_flex(v_idx), cot_rigid(v_idx), ...
            vel_flex_std(v_idx), vel_rigid_std(v_idx), ...
            accel_flex_std(v_idx), accel_rigid_std(v_idx), ...
            power_flex_std(v_idx), power_rigid_std(v_idx), ...
            cot_flex_std(v_idx), cot_rigid_std(v_idx), ...
            target_time, colors);
    end
end

%% Calculate Ratio Error Bars
% For ratio R = A/B: σ_R = R × sqrt((σ_A/A)² + (σ_B/B)²)
vel_ratio = vel_flex ./ vel_rigid;
vel_ratio_std = vel_ratio .* sqrt((vel_flex_std ./ vel_flex).^2 + (vel_rigid_std ./ vel_rigid).^2);

accel_ratio = accel_flex ./ accel_rigid;
accel_ratio_std = accel_ratio .* sqrt((accel_flex_std ./ accel_flex).^2 + (accel_rigid_std ./ accel_rigid).^2);

%% Figure 9: Ratio Plots (Voltage)
fig_ratio = figure('Position', [300, 300, 1200, 700], 'Color', 'w');
set(fig_ratio, 'DefaultAxesFontSize', 11, 'DefaultAxesFontName', 'Arial');

% Velocity ratio
subplot(2, 1, 1);
plot_ratio_with_error(voltages, vel_ratio, vel_ratio_std, colors.flex, ...
    'Velocity Ratio (Flex/Rigid)', ...
    sprintf('Velocity Ratio at t=%.1fs', target_time));

% Acceleration ratio
subplot(2, 1, 2);
plot_ratio_with_error(voltages, accel_ratio, accel_ratio_std, colors.rigid, ...
    'Acceleration Ratio (Flex/Rigid)', ...
    'Acceleration Ratio');

sgtitle('Performance Ratios: Flex 15mm / Rigid', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 10: Ratio Plots (Frequency)
fig_ratio_freq = figure('Position', [350, 350, 1200, 700], 'Color', 'w');
set(fig_ratio_freq, 'DefaultAxesFontSize', 11, 'DefaultAxesFontName', 'Arial');

% Velocity ratio vs frequency
subplot(2, 1, 1);
plot_ratio_with_error(frequencies, vel_ratio, vel_ratio_std, colors.flex, ...
    'Velocity Ratio (Flex/Rigid)', ...
    sprintf('Velocity Ratio at t=%.1fs', target_time), 'Frequency (Hz)');

% Acceleration ratio vs frequency
subplot(2, 1, 2);
plot_ratio_with_error(frequencies, accel_ratio, accel_ratio_std, colors.rigid, ...
    'Acceleration Ratio (Flex/Rigid)', ...
    'Acceleration Ratio', 'Frequency (Hz)');

sgtitle('Performance Ratios vs Frequency: Flex 15mm / Rigid', 'FontSize', 14, 'FontWeight', 'bold');

%% Print Results
print_results(fit_params, voltages, voltage_labels, vel_flex, vel_rigid, ...
    accel_flex, accel_rigid, power_flex, power_rigid, cot_flex, cot_rigid, target_time);

% Print fit quality summary
fprintf('\n========================================\n');
fprintf('FIT QUALITY (R² values)\n');
fprintf('========================================\n');
for v_idx = 1:length(voltages)
    flex_field = strrep(sprintf('Flex15mm_%.1fV', voltages(v_idx)), '.', '_');
    rigid_field = strrep(sprintf('Rigid_%.1fV', voltages(v_idx)), '.', '_');

    if isfield(fit_params, flex_field) && isfield(fit_params.(flex_field), 'r_squared')
        r2_flex = fit_params.(flex_field).r_squared;
        fprintf('%s Flex:  R² = %.5f', voltage_labels{v_idx}, r2_flex);
        if r2_flex < 0.99
            fprintf(' ⚠');
        end
        fprintf('\n');
    end

    if isfield(fit_params, rigid_field) && isfield(fit_params.(rigid_field), 'r_squared')
        r2_rigid = fit_params.(rigid_field).r_squared;
        fprintf('%s Rigid: R² = %.5f', voltage_labels{v_idx}, r2_rigid);
        if r2_rigid < 0.99
            fprintf(' ⚠');
        end
        fprintf('\n');
    end
end
fprintf('\nNote: R² > 0.99 indicates excellent fit\n');

fprintf('\n✓ Analysis complete! Generated 10 figures.\n');

%% ==================== HELPER FUNCTIONS ====================

function fit_params = process_condition(condition, voltage, files, ...
    traj_folder, sampling_freq, pixels_per_mm, bl_mm, target_time, poly_order, colors, fit_params)
    % Process trajectory data for one condition and fit polynomial curve with high accuracy

    % Format voltage for filename matching
    if voltage == 9.0
        voltage_str = '9V';  % Special case: 9V not 9.0V
    else
        voltage_str = sprintf('%.1fV', voltage);
    end
    pattern = sprintf('NoWire_%s_%s', condition, voltage_str);

    all_time = [];
    all_position = [];
    trial_velocities = [];  % Store velocity at target_time for each trial
    trial_accels = [];      % Store acceleration at target_time for each trial

    % Load all trials
    for trial = 1:3
        file_pattern = sprintf('%s_T%d', pattern, trial);
        file_idx = find(contains({files.name}, file_pattern), 1);

        if ~isempty(file_idx)
            data = readmatrix(fullfile(traj_folder, files(file_idx).name));
            X = data(~isnan(data(:,1)), 1);

            if ~isempty(X)
                % Process trajectory
                X = -X - min(-X);
                threshold_idx = find(X >= 30, 1);
                if ~isempty(threshold_idx)
                    X = movmedian(X(threshold_idx:end), 15) / pixels_per_mm;
                    time = (0:length(X)-1)' / sampling_freq;

                    % Fit individual trial to calculate per-trial metrics
                    if max(time) >= target_time
                        % Normalize individual trial
                        t_trial_mean = mean(time);
                        t_trial_std = std(time);
                        x_trial_mean = mean(X);
                        x_trial_std = std(X);

                        t_trial_norm = (time - t_trial_mean) / t_trial_std;
                        x_trial_norm = (X - x_trial_mean) / x_trial_std;

                        % Fit trial
                        p_trial_norm = polyfit(t_trial_norm, x_trial_norm, poly_order);

                        % Calculate velocity at target_time
                        t_target_norm_trial = (target_time - t_trial_mean) / t_trial_std;
                        p_vel_trial = polyder(p_trial_norm);
                        v_norm_trial = polyval(p_vel_trial, t_target_norm_trial);
                        v_trial = v_norm_trial * (x_trial_std / t_trial_std);
                        trial_velocities = [trial_velocities; v_trial / bl_mm];

                        % Calculate acceleration at target_time
                        if poly_order >= 2
                            p_accel_trial = polyder(p_vel_trial);
                            a_norm_trial = polyval(p_accel_trial, t_target_norm_trial);
                            a_trial = a_norm_trial * (x_trial_std / (t_trial_std^2));
                            trial_accels = [trial_accels; a_trial / bl_mm];
                        end
                    end

                    % Plot trial data
                    plot(time, X, '-', 'LineWidth', 1, 'Color', colors.gray);
                    hold on;

                    all_time = [all_time; time];
                    all_position = [all_position; X];
                end
            end
        end
    end

    % Fit polynomial curve with improved accuracy
    if ~isempty(all_time)
        [all_time, idx] = sort(all_time);
        all_position = all_position(idx);

        % Use polyfit with higher precision
        % For better numerical stability, center and scale the data
        t_mean = mean(all_time);
        t_std = std(all_time);
        x_mean = mean(all_position);
        x_std = std(all_position);

        % Normalized data for better conditioning
        t_norm = (all_time - t_mean) / t_std;
        x_norm = (all_position - x_mean) / x_std;

        % Fit polynomial with warning suppression for rank deficiency
        warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
        [p_norm, S, mu] = polyfit(t_norm, x_norm, poly_order);
        warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');

        % Transform back to original scale
        % For display, evaluate on normalized then transform back
        t_fit_norm = linspace(min(t_norm), max(t_norm), 200);
        x_fit_norm = polyval(p_norm, t_fit_norm);

        % Transform back
        t_fit = t_fit_norm * t_std + t_mean;
        x_fit = x_fit_norm * x_std + x_mean;

        % Calculate fit quality (R-squared)
        x_pred = polyval(p_norm, t_norm) * x_std + x_mean;
        ss_res = sum((all_position - x_pred).^2);
        ss_tot = sum((all_position - mean(all_position)).^2);
        r_squared = 1 - ss_res/ss_tot;

        % Plot fitted curve
        if contains(condition, 'Flex')
            plot(t_fit, x_fit, '-', 'LineWidth', 2.5, 'Color', colors.flex);
        else
            plot(t_fit, x_fit, '-', 'LineWidth', 2.5, 'Color', colors.rigid);
        end

        % Calculate derivatives in normalized space, then transform
        % For polynomial p(t), derivative coefficients are [n*p(1), (n-1)*p(2), ..., p(n-1)]
        field_name = strrep(sprintf('%s_%.1fV', condition, voltage), '.', '_');

        % Store R-squared for fit quality
        fit_params.(field_name).r_squared = r_squared;

        % CORRECTED: Use average of individual trial velocities, not combined fit
        if target_time <= max(all_time) && length(trial_velocities) >= 1
            % Method: Average individual trial velocities (CORRECT)
            v_at_target_bl = mean(trial_velocities);  % BL/s
            fit_params.(field_name).v_at_target_bl = v_at_target_bl;

            % Method: Average individual trial accelerations (CORRECT)
            if length(trial_accels) >= 1
                accel_bl = mean(trial_accels);  % BL/s²
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

            % Also calculate from combined fit for display (but don't use for comparison)
            t_target_norm = (target_time - t_mean) / t_std;
            p_vel_norm = polyder(p_norm);
            v_norm = polyval(p_vel_norm, t_target_norm);
            v_combined_fit = v_norm * (x_std / t_std) / bl_mm;

            % Initial velocity (at t=0)
            t_zero_norm = (0 - t_mean) / t_std;
            v_zero_norm = polyval(p_vel_norm, t_zero_norm);
            v_initial = v_zero_norm * (x_std / t_std);  % mm/s
            fit_params.(field_name).v_initial_bl = v_initial / bl_mm;  % BL/s

            % Add annotation with fit quality and std dev
            if length(trial_velocities) >= 2
                text(0.98, 0.98, sprintf('v: %.2f±%.2f BL/s\na: %.3f±%.3f BL/s²\nR²: %.4f', ...
                    v_at_target_bl, fit_params.(field_name).v_at_target_std_bl, ...
                    accel_bl, fit_params.(field_name).accel_std_bl, r_squared), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', 'BackgroundColor', 'w', ...
                    'EdgeColor', 'k', 'FontSize', 9, 'Margin', 3);
            else
                text(0.98, 0.98, sprintf('v: %.2f BL/s\na: %.3f BL/s²\nR²: %.4f', ...
                    v_at_target_bl, accel_bl, r_squared), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', 'BackgroundColor', 'w', ...
                    'EdgeColor', 'k', 'FontSize', 9, 'Margin', 3);
            end
        end
    end

    % Format subplot
    xlabel('Time (s)', 'FontWeight', 'bold');
    ylabel('Position (mm)', 'FontWeight', 'bold');
    title(strrep(sprintf('%s %.1fV', condition, voltage), '_', ' '), 'FontSize', 11);
    grid off; box on;
    if ~isempty(all_time)
        xlim([0 max(all_time)]);
    end
end

function [vel_flex, vel_rigid, accel_flex, accel_rigid, vel_flex_std, vel_rigid_std, accel_flex_std, accel_rigid_std] = extract_metrics(fit_params, voltages)
    % Extract velocity and acceleration arrays from fit parameters with std dev

    vel_flex = nan(size(voltages));
    vel_rigid = nan(size(voltages));
    accel_flex = nan(size(voltages));
    accel_rigid = nan(size(voltages));
    vel_flex_std = nan(size(voltages));
    vel_rigid_std = nan(size(voltages));
    accel_flex_std = nan(size(voltages));
    accel_rigid_std = nan(size(voltages));

    for i = 1:length(voltages)
        flex_field = strrep(sprintf('Flex15mm_%.1fV', voltages(i)), '.', '_');
        rigid_field = strrep(sprintf('Rigid_%.1fV', voltages(i)), '.', '_');

        if isfield(fit_params, flex_field)
            if isfield(fit_params.(flex_field), 'v_at_target_bl')
                vel_flex(i) = fit_params.(flex_field).v_at_target_bl;
            end
            if isfield(fit_params.(flex_field), 'accel_bl')
                accel_flex(i) = fit_params.(flex_field).accel_bl;
            end
            if isfield(fit_params.(flex_field), 'v_at_target_std_bl')
                vel_flex_std(i) = fit_params.(flex_field).v_at_target_std_bl;
            end
            if isfield(fit_params.(flex_field), 'accel_std_bl')
                accel_flex_std(i) = fit_params.(flex_field).accel_std_bl;
            end
        end

        if isfield(fit_params, rigid_field)
            if isfield(fit_params.(rigid_field), 'v_at_target_bl')
                vel_rigid(i) = fit_params.(rigid_field).v_at_target_bl;
            end
            if isfield(fit_params.(rigid_field), 'accel_bl')
                accel_rigid(i) = fit_params.(rigid_field).accel_bl;
            end
            if isfield(fit_params.(rigid_field), 'v_at_target_std_bl')
                vel_rigid_std(i) = fit_params.(rigid_field).v_at_target_std_bl;
            end
            if isfield(fit_params.(rigid_field), 'accel_std_bl')
                accel_rigid_std(i) = fit_params.(rigid_field).accel_std_bl;
            end
        end
    end
end

function plot_comparison_with_error(voltages, data_flex, data_rigid, std_flex, std_rigid, colors, ylabel_str, title_str)
    % Create comparison line plot with error bars

    % Plot with error bars
    errorbar(voltages, data_flex, std_flex, 'o-', 'LineWidth', 2.5, 'MarkerSize', 9, ...
        'Color', colors.flex, 'MarkerFaceColor', colors.flex, 'DisplayName', 'Flex 15mm', ...
        'CapSize', 8, 'LineStyle', '-');
    hold on;
    errorbar(voltages, data_rigid, std_rigid, 's-', 'LineWidth', 2.5, 'MarkerSize', 9, ...
        'Color', colors.rigid, 'MarkerFaceColor', colors.rigid, 'DisplayName', 'Rigid', ...
        'CapSize', 8, 'LineStyle', '-');

    xlabel('Voltage (V)', 'FontWeight', 'bold');
    ylabel(ylabel_str, 'FontWeight', 'bold');
    title(title_str, 'FontSize', 12);
    legend('Location', 'northwest', 'FontSize', 10);
    grid off; box on;
    xlim([min(voltages)-0.2, max(voltages)+0.2]);
end

function create_voltage_figure(v_idx, voltage_label, vel_flex, vel_rigid, ...
    accel_flex, accel_rigid, power_flex, power_rigid, cot_flex, cot_rigid, ...
    vel_flex_std, vel_rigid_std, accel_flex_std, accel_rigid_std, ...
    power_flex_std, power_rigid_std, cot_flex_std, cot_rigid_std, ...
    target_time, colors)
    % Create individual voltage comparison figure with error bars (4 panels)

    fig = figure('Position', [200 + v_idx*30, 200 + v_idx*30, 900, 800], 'Color', 'w');
    set(fig, 'DefaultAxesFontSize', 11, 'DefaultAxesFontName', 'Arial');

    % Velocity bar chart
    subplot(2, 2, 1);
    create_bar_subplot_with_error([vel_rigid, vel_flex], [vel_rigid_std, vel_flex_std], colors, ...
        sprintf('Velocity at t=%.1fs (BL/s)', target_time), ...
        sprintf('Velocity at %s', voltage_label));

    % Acceleration bar chart
    subplot(2, 2, 2);
    create_bar_subplot_with_error([accel_rigid, accel_flex], [accel_rigid_std, accel_flex_std], colors, ...
        'Acceleration (BL/s²)', ...
        sprintf('Acceleration at %s', voltage_label));

    % Power bar chart
    subplot(2, 2, 3);
    create_bar_subplot_with_error([power_rigid, power_flex], [power_rigid_std, power_flex_std], colors, ...
        'Power (mW)', ...
        sprintf('Power at %s', voltage_label));

    % Cost of Transport bar chart
    subplot(2, 2, 4);
    create_bar_subplot_with_error([cot_rigid, cot_flex], [cot_rigid_std, cot_flex_std], colors, ...
        'Cost of Transport (J/(N·m))', ...
        sprintf('Cost of Transport at %s', voltage_label));

    sgtitle(sprintf('%s: Flex 15mm vs Rigid', voltage_label), ...
        'FontSize', 13, 'FontWeight', 'bold');
end

function create_bar_subplot_with_error(data, std_data, colors, ylabel_str, title_str)
    % Create bar chart subplot with error bars and percentage annotation

    b = bar(data);
    b.FaceColor = 'flat';
    b.CData(1,:) = colors.rigid;
    b.CData(2,:) = colors.flex;
    b.EdgeColor = 'k';
    b.LineWidth = 1.2;

    % Add error bars
    hold on;
    x_pos = 1:length(data);
    errorbar(x_pos, data, std_data, 'k.', 'LineWidth', 1.5, 'CapSize', 10, ...
        'HandleVisibility', 'off');

    set(gca, 'XTickLabel', {'Rigid', 'Flex 15mm'}, 'FontSize', 11);
    ylabel(ylabel_str, 'FontWeight', 'bold');
    title(title_str, 'FontSize', 11);
    grid off; box on;

    % Add percentage
    pct = ((data(2) - data(1)) / data(1)) * 100;
    color_text = colors.green;
    if pct < 0
        color_text = colors.red;
    end
    ylims = ylim;
    text(1.5, ylims(2)*0.92, sprintf('%+.1f%%', pct), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
        'FontSize', 14, 'Color', color_text);
end

function plot_ratio_with_error(x_values, ratio, ratio_std, color, ylabel_str, title_str, xlabel_str)
    % Create ratio plot with baseline and error bars

    % Default xlabel
    if nargin < 7
        xlabel_str = 'Voltage (V)';
    end

    % Plot with error bars
    errorbar(x_values, ratio, ratio_std, 'o-', 'LineWidth', 2.5, 'MarkerSize', 10, ...
        'Color', color, 'MarkerFaceColor', color, 'CapSize', 8);
    hold on;
    yline(1, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Add percentage labels
    for i = 1:length(x_values)
        if ~isnan(ratio(i))
            pct = (ratio(i) - 1) * 100;
            text(x_values(i), ratio(i) + 0.015, sprintf('%+.1f%%', pct), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
        end
    end

    xlabel(xlabel_str, 'FontWeight', 'bold');
    ylabel(ylabel_str, 'FontWeight', 'bold');
    title(title_str, 'FontSize', 12);
    grid off; box on;
    xlim([min(x_values)-0.5, max(x_values)+0.5]);
end

function print_results(fit_params, voltages, voltage_labels, vel_flex, ...
    vel_rigid, accel_flex, accel_rigid, power_flex, power_rigid, cot_flex, cot_rigid, target_time)
    % Print analysis results to console

    fprintf('\n========================================\n');
    fprintf('PERFORMANCE COMPARISON SUMMARY\n');
    fprintf('========================================\n');
    fprintf('Analysis time point: t = %.1fs\n\n', target_time);

    fprintf('%-8s | %-10s | %-10s | %-10s | %-10s | %-10s\n', ...
        'Voltage', 'Vel Δ', 'Accel Δ', 'Power Δ', 'CoT Δ', 'Status');
    fprintf('---------|------------|------------|------------|------------|------------\n');

    pct_vel = ((vel_flex - vel_rigid) ./ vel_rigid) * 100;
    pct_accel = ((accel_flex - accel_rigid) ./ accel_rigid) * 100;
    pct_power = ((power_flex - power_rigid) ./ power_rigid) * 100;
    pct_cot = ((cot_flex - cot_rigid) ./ cot_rigid) * 100;

    for i = 1:length(voltages)
        status = '✓';
        if pct_vel(i) < 0 || pct_accel(i) < 0
            status = '⚠';
        end
        fprintf('%-8s | %+9.1f%% | %+9.1f%% | %+9.1f%% | %+9.1f%% | %s\n', ...
            voltage_labels{i}, pct_vel(i), pct_accel(i), pct_power(i), pct_cot(i), status);
    end

    fprintf('\n');
    fprintf('Average improvements:\n');
    fprintf('  Velocity:     %+.1f%%\n', mean(pct_vel, 'omitnan'));
    fprintf('  Acceleration: %+.1f%%\n', mean(pct_accel, 'omitnan'));
    fprintf('  Power:        %+.1f%%\n', mean(pct_power, 'omitnan'));
    fprintf('  CoT:          %+.1f%%\n', mean(pct_cot, 'omitnan'));
    fprintf('\n');
end

function [power_flex, power_rigid, power_flex_std, power_rigid_std] = load_power_data(power_folder, voltages, conditions)
    % Load and average power consumption data for each condition
    % Reads CROPPED power data (CSV format: time, power)
    % Returns power in mW

    power_flex = nan(size(voltages));
    power_rigid = nan(size(voltages));
    power_flex_std = nan(size(voltages));
    power_rigid_std = nan(size(voltages));

    for c_idx = 1:length(conditions)
        condition = conditions{c_idx};

        for v_idx = 1:length(voltages)
            voltage = voltages(v_idx);

            % Format voltage for filename matching
            if voltage == 9.0
                voltage_str = '9V';
            else
                voltage_str = sprintf('%.1fV', voltage);
            end

            % Load all 3 trials
            trial_powers = [];

            for trial = 1:3
                % Construct filename for CROPPED data
                % Format: Cropped_[Nozzle]_[Voltage]_T[Trial].txt
                filename = sprintf('Cropped_%s_%s_T%d.txt', condition, strrep(voltage_str, '.', '_'), trial);
                filepath = fullfile(power_folder, filename);

                if isfile(filepath)
                    try
                        % Read cropped CSV file (columns: time, power)
                        data = readmatrix(filepath);

                        if ~isempty(data) && size(data, 2) >= 2
                            time = data(:, 1);
                            power = data(:, 2);

                            % Average power over entire cropped period
                            % (cropping already removed startup/shutdown transients)
                            trial_powers(end+1) = mean(power);
                        end
                    catch ME
                        fprintf('Warning: Error reading %s: %s\n', filename, ME.message);
                    end
                end
            end

            % Calculate mean and std across trials
            if ~isempty(trial_powers)
                if strcmp(condition, 'Flex15mm')
                    power_flex(v_idx) = mean(trial_powers);
                    if length(trial_powers) >= 2
                        power_flex_std(v_idx) = std(trial_powers);
                    else
                        power_flex_std(v_idx) = 0;
                    end
                else
                    power_rigid(v_idx) = mean(trial_powers);
                    if length(trial_powers) >= 2
                        power_rigid_std(v_idx) = std(trial_powers);
                    else
                        power_rigid_std(v_idx) = 0;
                    end
                end
            end
        end
    end
end
