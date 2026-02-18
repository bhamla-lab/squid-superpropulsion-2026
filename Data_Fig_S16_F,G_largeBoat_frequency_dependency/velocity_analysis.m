% Velocity analysis using polynomial fitting
% Plots fitted velocity curves for all conditions

% METHOD:
%   1. Normalizes data (center & scale) for numerical stability
%   2. Fits polynomial of adjustable order (default: cubic)
%   3. Calculates velocity analytically using polyder()
%   4. Transforms results back to original scale
%   5. Plots velocity curves for each trial and combined fit

clear; close all; clc;

% PARAMETERS
target_time = 5;        % Time point in data for velocity/acceleration analysis
sampling_freq = 30;     % Hz/video frame rate
pixels_per_mm = 2.9;    % Calibration from DLTdv8
bl_mm = 260;            % Boat length in mm
poly_order = 3;         % Polynomial order for fitting

% Folder info
traj_folder = 'Trajectories';
voltages = [7.8, 8.1, 8.4, 8.7, 9.0, 9.3];
frequencies = [16, 17, 18, 21, 22, 24];  % Hz corresponding to each voltage
voltage_labels = arrayfun(@(v) sprintf('%.1fV', v), voltages, 'UniformOutput', false);
conditions = {'Flex15mm', 'Rigid'};

% Color scheme
colors.flex = [0.2 0.6 0.8];
colors.rigid = [0.8 0.4 0.2];
colors.gray = [0.7 0.7 0.7];

%% Load and Process Data
fprintf('Processing trajectory data for velocity plots...\n');

% Get trajectory files
files = dir(fullfile(traj_folder, '*xypts.csv'));
files = files(~startsWith({files.name}, '._'));

%% Figure 1: Velocity Data with Fitted Curves
fig1 = figure('Position', [100, 100, 1800, 1200], 'Color', 'w');
set(fig1, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Arial');

plot_idx = 1;
for cond = conditions
    for v_idx = 1:length(voltages)
        subplot(3, 4, plot_idx);
        process_velocity_condition(cond{1}, voltages(v_idx), files, ...
            traj_folder, sampling_freq, pixels_per_mm, bl_mm, target_time, ...
            poly_order, colors);
        plot_idx = plot_idx + 1;
    end
end

sgtitle('NoWire Validation - Fitted Velocity Curves', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\n Analysis complete! Generated velocity figure.\n');

%% ==================== HELPER FUNCTIONS ====================

function process_velocity_condition(condition, voltage, files, ...
    traj_folder, sampling_freq, pixels_per_mm, bl_mm, target_time, poly_order, colors)
    % Process trajectory data for one condition and plot fitted velocity curves

    % Format voltage for filename matching
    if voltage == 9.0
        voltage_str = '9V';  % Special case: 9V not 9.0V
    else
        voltage_str = sprintf('%.1fV', voltage);
    end
    pattern = sprintf('NoWire_%s_%s', condition, voltage_str);

    all_time = [];
    all_position = [];
    trial_data = {};  % Store individual trial data for velocity plotting

    trial_velocities_at_target = [];

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

                    % Store trial data
                    trial_data{end+1} = struct('time', time, 'position', X);

                    % Fit individual trial and calculate velocity curve
                    t_trial_mean = mean(time);
                    t_trial_std = std(time);
                    x_trial_mean = mean(X);
                    x_trial_std = std(X);

                    t_trial_norm = (time - t_trial_mean) / t_trial_std;
                    x_trial_norm = (X - x_trial_mean) / x_trial_std;

                    % Fit polynomial
                    warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
                    p_trial_norm = polyfit(t_trial_norm, x_trial_norm, poly_order);
                    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');

                    % Calculate velocity curve using polyder
                    p_vel_norm = polyder(p_trial_norm);

                    % Evaluate velocity over time range
                    t_fit_norm = linspace(min(t_trial_norm), max(t_trial_norm), 200);
                    v_fit_norm = polyval(p_vel_norm, t_fit_norm);

                    % Transform back to original scale
                    t_fit = t_fit_norm * t_trial_std + t_trial_mean;
                    v_fit = v_fit_norm * (x_trial_std / t_trial_std);  % mm/s
                    v_fit_bl = v_fit / bl_mm;  % BL/s

                    % Plot individual trial velocity
                    plot(t_fit, v_fit_bl, '-', 'LineWidth', 1, 'Color', colors.gray);
                    hold on;

                    % Store velocity at target time
                    if max(time) >= target_time
                        t_target_norm = (target_time - t_trial_mean) / t_trial_std;
                        v_at_target = polyval(p_vel_norm, t_target_norm) * (x_trial_std / t_trial_std) / bl_mm;
                        trial_velocities_at_target = [trial_velocities_at_target; v_at_target];
                    end

                    all_time = [all_time; time];
                    all_position = [all_position; X];
                end
            end
        end
    end

    % Fit combined polynomial and plot velocity curve
    if ~isempty(all_time)
        [all_time_sorted, idx] = sort(all_time);
        all_position_sorted = all_position(idx);

        % Normalize combined data
        t_mean = mean(all_time_sorted);
        t_std = std(all_time_sorted);
        x_mean = mean(all_position_sorted);
        x_std = std(all_position_sorted);

        t_norm = (all_time_sorted - t_mean) / t_std;
        x_norm = (all_position_sorted - x_mean) / x_std;

        % Fit polynomial
        warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
        p_norm = polyfit(t_norm, x_norm, poly_order);
        warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');

        % Calculate velocity using polyder
        p_vel_norm = polyder(p_norm);

        % Evaluate over time range
        t_fit_norm = linspace(min(t_norm), max(t_norm), 200);
        v_fit_norm = polyval(p_vel_norm, t_fit_norm);

        % Transform back
        t_fit = t_fit_norm * t_std + t_mean;
        v_fit = v_fit_norm * (x_std / t_std);  % mm/s
        v_fit_bl = v_fit / bl_mm;  % BL/s

        % Plot fitted velocity curve (combined)
        if contains(condition, 'Flex')
            plot(t_fit, v_fit_bl, '-', 'LineWidth', 2.5, 'Color', colors.flex);
        else
            plot(t_fit, v_fit_bl, '-', 'LineWidth', 2.5, 'Color', colors.rigid);
        end

        % Calculate R-squared for position fit
        x_pred = polyval(p_norm, t_norm) * x_std + x_mean;
        ss_res = sum((all_position_sorted - x_pred).^2);
        ss_tot = sum((all_position_sorted - mean(all_position_sorted)).^2);
        r_squared = 1 - ss_res/ss_tot;

        % Add vertical line at target time
        xline(target_time, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');

        % Calculate mean velocity at target time
        if ~isempty(trial_velocities_at_target)
            v_mean = mean(trial_velocities_at_target);
            v_std = std(trial_velocities_at_target);

            if length(trial_velocities_at_target) >= 2
                text(0.98, 0.98, sprintf('v@%.0fs: %.2f +/- %.2f BL/s\nR^2: %.4f', ...
                    target_time, v_mean, v_std, r_squared), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', 'BackgroundColor', 'w', ...
                    'EdgeColor', 'k', 'FontSize', 9, 'Margin', 3);
            else
                text(0.98, 0.98, sprintf('v@%.0fs: %.2f BL/s\nR^2: %.4f', ...
                    target_time, v_mean, r_squared), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', 'BackgroundColor', 'w', ...
                    'EdgeColor', 'k', 'FontSize', 9, 'Margin', 3);
            end
        end
    end

    % Format subplot
    xlabel('Time (s)', 'FontWeight', 'bold');
    ylabel('Velocity (BL/s)', 'FontWeight', 'bold');
    title(strrep(sprintf('%s %.1fV', condition, voltage), '_', ' '), 'FontSize', 11);
    grid off; box on;
    if ~isempty(all_time)
        xlim([0 max(all_time)]);
        ylim([0 inf]);  % Velocity should be positive
    end
end
