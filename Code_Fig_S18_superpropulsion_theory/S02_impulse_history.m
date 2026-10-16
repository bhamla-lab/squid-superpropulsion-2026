%% S02_impulse_history.m
% Compare predicted impulse history with experimental data
% Case: 30WT 0200MS
% Models: Base and Collapse
%
% Theory:
%   Hydrodynamic impulse is proportional to the cumulative integral of u^2:
%     I_h(t) ~ integral_0^t u(t')^2 dt'
%   Each model's impulse is normalized by the total input impulse I_in.
%   Experimental impulse is normalized by the 00WT plateau value.

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

% Input pulse duration T [s] for 0200MS
T_input = 0.106532;

%% ========================================================================
%  Derived quantities
%  ========================================================================
c_wave     = sqrt(Eh_30WT / (2 * rho * R));
tau        = L / c_wave;
tau_T      = tau / T_input;
omega0_norm = pi / (2 * tau_T);

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

% Total injected volume (normalized)
mask_in    = V_norm > 0;
Q_in_total = sum(V_norm(mask_in)) * dt_norm;

% Cumulative input impulse: I_in ~ integral(u^2 dt)
I_in_cum   = cumsum(V_norm.^2) * dt_norm;
I_in_total = I_in_cum(end);

% Experimental impulse data (trimmed)
imp = load(fullfile(data_dir, 'impulse_exp.mat'));
t_input_raw = imp.t_input;
I_input_raw = imp.I_h_input;
t_exp       = imp.t_exp;
I_h_exp_raw = imp.I_h_exp;
I_h_00WT    = imp.I_h_00WT;   % normalization constant

fprintf('I_h_00WT (0200MS) = %.6f kg*m/s\n', I_h_00WT);

%% ========================================================================
%  Solve ODE: q_out'' = omega0^2 * (V_in - q_out)
%  ========================================================================
q_out     = zeros(n, 1);
q_out_dot = 0;

for i = 2:n
    q_out_ddot = omega0_norm^2 * (V_norm(i-1) - q_out(i-1));
    q_out_dot  = q_out_dot + q_out_ddot * dt_norm;
    q_out(i)   = q_out(i-1) + q_out_dot * dt_norm;
end

%% ========================================================================
%  Base Model
%  ========================================================================
first_neg_idx = find(q_out < 0, 1, 'first');
if isempty(first_neg_idx), first_neg_idx = n + 1; end

q_base = q_out;
q_base(first_neg_idx:end) = 0;

%% ========================================================================
%  Collapse Model
%  ========================================================================
Q_out_cumsum  = cumsum(q_out) * dt_norm;
vol_match_idx = find(Q_out_cumsum >= Q_in_total, 1, 'first');
if isempty(vol_match_idx), vol_match_idx = n; end

q_0 = q_out(vol_match_idx);

q_collapse = zeros(n, 1);
if q_0 > 0
    Del_t_vol   = Q_in_total / q_0;
    n_plat_vol  = round(Del_t_vol / dt_norm);
    end_vol_idx = min(vol_match_idx + n_plat_vol, n);
    q_collapse(1:vol_match_idx)             = q_out(1:vol_match_idx);
    q_collapse(vol_match_idx+1:end_vol_idx) = q_0;
end

%% ========================================================================
%  Cumulative impulse: I ~ integral(u^2 dt), normalized by I_in_total
%  ========================================================================
I_in_norm       = I_in_cum / I_in_total;
I_base_cum      = cumsum(q_base.^2) * dt_norm / I_in_total;
I_collapse_cum  = cumsum(q_collapse.^2) * dt_norm / I_in_total;

%% ========================================================================
%  Plot
%  ========================================================================
scatter_interval = 1;

figure('Position', [200 200 600 450]);
hold on;

% Theory curves
h1 = plot(t_norm, I_in_norm, 'k-', 'LineWidth', 1.5);
h2 = plot(t_norm, I_base_cum, 'r--', 'LineWidth', 1.5);
h3 = plot(t_norm, I_collapse_cum, 'b-', 'LineWidth', 1.5);

% Experimental scatter
idx_in  = 1:scatter_interval:length(t_input_raw);
idx_exp = 1:scatter_interval:length(t_exp);
h4 = scatter(t_input_raw(idx_in), I_input_raw(idx_in) / I_h_00WT, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
h5 = scatter(t_exp(idx_exp), I_h_exp_raw(idx_exp) / I_h_00WT, 20, 'r', 'filled');

xlabel('$\textit{\textsf{t/T}}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\textit{\textsf{I}}_{\textit{\textsf{h}}} / \textit{\textsf{I}}_{\textit{\textsf{h,00WT}}}$', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('Impulse History: 30WT 0200MS (\\tau/T = %.2f)', tau_T));
xlim([0 3]);
ylim([-0.5 5]);
legend([h1 h2 h3 h4 h5], {'I_{in} (fit)', 'Base', 'Collapse', ...
    'I_{in} (rigid)', 'Exp (flexible)'}, 'Location', 'southeast', 'FontSize', 8);

