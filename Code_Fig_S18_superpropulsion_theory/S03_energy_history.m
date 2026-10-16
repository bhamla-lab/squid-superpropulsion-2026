%% S03_energy_history.m
% Compare predicted kinetic energy history with experimental data
% Case: 30WT 0200MS
% Models: Base and Collapse
%
% Theory:
%   Kinetic energy flux ~ 1/2 * rho * u^2 * u * A, so cumulative KE ~ integral(u^3 dt)
%   Each model's KE is normalized by the total input KE.
%   Experimental KE is normalized by the 00WT plateau value.

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

% Cumulative input KE: KE_in ~ integral(u^3 dt)
KE_in_cum   = cumsum(V_norm.^3) * dt_norm;
KE_in_total = KE_in_cum(end);

% Experimental energy data (trimmed)
ke = load(fullfile(data_dir, 'energy_exp.mat'));
t_input_raw  = ke.t_input;
KE_input_raw = ke.KE_input;
t_exp        = ke.t_exp;
KE_exp_raw   = ke.KE_exp;
KE_00WT      = ke.KE_00WT;   % normalization constant

fprintf('KE_00WT (0200MS) = %.6f J\n', KE_00WT);

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
%  Cumulative KE: KE ~ integral(u^3 dt), normalized by KE_in_total
%  ========================================================================
KE_in_norm       = KE_in_cum / KE_in_total;
KE_base_cum      = cumsum(q_base.^3) * dt_norm / KE_in_total;
KE_collapse_cum  = cumsum(q_collapse.^3) * dt_norm / KE_in_total;

%% ========================================================================
%  Plot
%  ========================================================================
scatter_interval = 1;

figure('Position', [200 200 600 450]);
hold on;

% Theory curves
h1 = plot(t_norm, KE_in_norm, 'k-', 'LineWidth', 1.5);
h2 = plot(t_norm, KE_base_cum, 'r--', 'LineWidth', 1.5);
h3 = plot(t_norm, KE_collapse_cum, 'b-', 'LineWidth', 1.5);

% Experimental scatter
idx_in  = 1:scatter_interval:length(t_input_raw);
idx_exp = 1:scatter_interval:length(t_exp);
h4 = scatter(t_input_raw(idx_in), KE_input_raw(idx_in) / KE_00WT, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
h5 = scatter(t_exp(idx_exp), KE_exp_raw(idx_exp) / KE_00WT, 20, 'r', 'filled');

xlabel('$\textit{\textsf{t/T}}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\textit{\textsf{KE}} / \textit{\textsf{KE}}_{\textit{\textsf{00WT}}}$', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('Kinetic Energy History: 30WT 0200MS (\\tau/T = %.2f)', tau_T));
xlim([0 3]);
ylim([-0.5 10]);
legend([h1 h2 h3 h4 h5], {'KE_{in} (fit)', 'Base', 'Collapse', ...
    'KE_{in} (rigid)', 'Exp (flexible)'}, 'Location', 'southeast', 'FontSize', 8);

