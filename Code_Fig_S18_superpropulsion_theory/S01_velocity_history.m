%% S01_velocity_history.m
% Compare predicted velocity history with experimental data
% Case: 30WT 0200MS
% Models: Base and Collapse
%
% Theory:
%   The nozzle exit velocity q_out is governed by a 2nd-order ODE:
%     q_out'' = omega0^2 * (V_in - q_out)
%   where omega0 = pi / (2*tau), tau = L / c_wave, c_wave = sqrt(Eh/(2*rho*R))
%
%   Base model:     raw ODE output, clipped to zero at first negative crossing
%   Collapse model: ODE output up to volume-match point, then constant plateau

clear; clc; close all;

%% ========================================================================
%  Physical Parameters
%  ========================================================================
D   = 15e-3;       % nozzle diameter [m]
R   = D / 2;       % nozzle radius [m]
L   = 2 * D;       % tube length [m]
rho = 1000;        % fluid density [kg/m^3]

% Wall stiffness Eh [N/m] for 30WT
Eh_30WT = 14.4;

% Input pulse duration T [s] for 0200MS (from C10a normalization)
T_input = 0.106532;

%% ========================================================================
%  Derived quantities
%  ========================================================================
% Wave speed in compliant tube
c_wave = sqrt(Eh_30WT / (2 * rho * R));

% Delay time tau [s]
tau = L / c_wave;

% Normalized delay
tau_T = tau / T_input;

% Natural frequency (normalized)
omega0_norm = pi / (2 * tau_T);

fprintf('30WT 0200MS parameters:\n');
fprintf('  c_wave  = %.3f m/s\n', c_wave);
fprintf('  tau     = %.4f s\n', tau);
fprintf('  tau/T   = %.3f\n', tau_T);
fprintf('  omega0  = %.3f (normalized)\n', omega0_norm);

%% ========================================================================
%  Load data (local, self-contained)
%  ========================================================================
data_dir = fullfile(fileparts(mfilename('fullpath')), 'data');

% Fitted input curve
inp = load(fullfile(data_dir, 'input_curve.mat'));
t_norm = inp.t_norm;
V_norm = inp.V_norm;

dt_norm = t_norm(2) - t_norm(1);
n       = length(t_norm);

% Total injected volume (normalized) - used for volume matching
mask_in    = V_norm > 0;
Q_in_total = sum(V_norm(mask_in)) * dt_norm;

% Experimental velocity data (already normalized and trimmed)
vel = load(fullfile(data_dir, 'velocity_exp.mat'));
t_input_raw = vel.t_input;
V_input_raw = vel.V_input;
t_exp       = vel.t_exp;
V_exp       = vel.V_exp;

%% ========================================================================
%  Solve ODE: q_out'' = omega0^2 * (V_in - q_out)
%  Forward Euler time integration
%  ========================================================================
q_out     = zeros(n, 1);
q_out_dot = 0;               % initial velocity of q_out

for i = 2:n
    q_out_ddot = omega0_norm^2 * (V_norm(i-1) - q_out(i-1));
    q_out_dot  = q_out_dot + q_out_ddot * dt_norm;
    q_out(i)   = q_out(i-1) + q_out_dot * dt_norm;
end

%% ========================================================================
%  Base Model: clip at first negative crossing
%  ========================================================================
first_neg_idx = find(q_out < 0, 1, 'first');
if isempty(first_neg_idx), first_neg_idx = n + 1; end

q_base = q_out;
q_base(first_neg_idx:end) = 0;

%% ========================================================================
%  Collapse Model: plateau at volume-matching point
%  ========================================================================
Q_out_cumsum  = cumsum(q_out) * dt_norm;
vol_match_idx = find(Q_out_cumsum >= Q_in_total, 1, 'first');
if isempty(vol_match_idx), vol_match_idx = n; end

q_0 = q_out(vol_match_idx);   % velocity at volume-match point

q_collapse = zeros(n, 1);
if q_0 > 0
    Del_t_vol  = Q_in_total / q_0;
    n_plat_vol = round(Del_t_vol / dt_norm);
    end_vol_idx = min(vol_match_idx + n_plat_vol, n);
    q_collapse(1:vol_match_idx)         = q_out(1:vol_match_idx);
    q_collapse(vol_match_idx+1:end_vol_idx) = q_0;
end

%% ========================================================================
%  Plot
%  ========================================================================
scatter_interval = 3;

figure('Position', [200 200 600 450]);
hold on;

% Theory curves
h1 = plot(t_norm, V_norm, 'k-', 'LineWidth', 1.5);
h2 = plot(t_norm, q_base, 'r--', 'LineWidth', 1.5);
h3 = plot(t_norm, q_collapse, 'b-', 'LineWidth', 1.5);

% Experimental scatter
idx_in  = 1:scatter_interval:length(t_input_raw);
idx_exp = 1:scatter_interval:length(t_exp);
h4 = scatter(t_input_raw(idx_in), V_input_raw(idx_in), 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
h5 = scatter(t_exp(idx_exp), V_exp(idx_exp), 20, 'r', 'filled');

plot([0 3], [0 0], 'k:', 'LineWidth', 0.5);

xlabel('$\textit{\textsf{t/T}}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\textit{\textsf{u/u}}_{\textit{\textsf{max}}}$', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('Velocity History: 30WT 0200MS (\\tau/T = %.2f)', tau_T));
xlim([0 1.5]);
ylim([-0.5 2.5]);
legend([h1 h2 h3 h4 h5], {'Input (fit)', 'Base', 'Collapse', ...
    'Input (rigid)', 'Exp (flexible)'}, 'Location', 'northeast', 'FontSize', 8);

