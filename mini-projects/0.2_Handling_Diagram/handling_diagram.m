% Handling Diagram
% At steady-state, set derivatives to zero. 
% Solve the linear system algebraically for yaw rate and sideslip vs speed and steering. 
% Gives the yaw rate gain plot — the most fundamental characterisation of vehicle handling.
%
% In the section 9 of the linear_bicycle_model.pdf file, about state-space form, we have the state vector [v ; r]
% This state v is set to 0 for the steadystate turning behavior Of the 2-DOF vehicle.
%
%
% [v_ss; r_ss] = -A \ B * delta        (per unit delta)
%
% Outputs per unit steer (delta = 1 rad):
%    r/delta         yaw rate gain           [rad/s per rad]
%    (1/R)/delta     curvature response      [1/m per rad]
%    (V^2/R)/delta   lateral accel response  [m/s^2 per rad]
%    beta/delta      sideslip response       [rad per rad]

%thisDir = fileparts(mfilename('fullpath'));
%addpath(fullfile(thisDir, '..', '..', '0.1_linear_bicycle_model', 'linear_bicycle_model.m'));

clear; clc; close all;

% ---- Base parameters ----
m  = 1500;   Iz = 2500;
a  = 1.2;    b  = 1.4;
L  = a + b;

% ---- Three vehicles ----
CFs = [100000,  80000*(b/a),   60000];
CRs = [ 60000,  80000,      100000];
labels = {'US','NS','OS'};
colors = {[0.20 0.60 0.90], [0.40 0.75 0.45], [0.90 0.40 0.35]};

K = m / L^2 * (a./CRs - b./CFs);

% ---- Sweep ----
V_vec = linspace(1, 60, 400);
nV = numel(V_vec);

G_r    = zeros(3, nV);
G_curv = zeros(3, nV);
G_ay   = zeros(3, nV);
G_beta = zeros(3, nV);

for i = 1:3
    for k = 1:nV
        V = V_vec(k);
        [A, B] = linear_bicycle_model(V, CFs(i), CRs(i));
        X = -A \ B;
        v_ss = X(1);
        r_ss = X(2);
        G_r(i,k)    = r_ss;
        G_curv(i,k) = r_ss / V;
        G_ay(i,k)   = V * r_ss;
        G_beta(i,k) = v_ss / V;
    end
end

% ---- Characteristic / critical speeds per vehicle ----
K_tol  = 1e-10;
V_ch   = nan(1,3);   % characteristic speed (only for K > 0)
V_cr   = nan(1,3);   % critical speed       (only for K < 0)

for i = 1:3
    if K(i) > K_tol
        V_ch(i) = sqrt(1 / K(i));
    elseif K(i) < -K_tol
        V_cr(i) = sqrt(-1 / K(i));
    end
    % K(i) ≈ 0 : neither exists, both stay NaN
end

% =====================================================================
% PLOT 1 : YAW RATE GAIN  r/delta
% =====================================================================
figure('Color','w','Name','Yaw rate gain');
hold on; grid on; box on;
for i = 1:3
    plot(V_vec, G_r(i,:), 'LineWidth', 1.8, 'Color', colors{i}, ...
         'DisplayName', sprintf('%s   K = %+.4f', labels{i}, K(i)));
end
for i = 1:3
    if ~isnan(V_ch(i))
        xline(V_ch(i), '--', 'Color', colors{i}, ...
              'Label', sprintf('V_{ch} = %.1f m/s', V_ch(i)), ...
              'LabelVerticalAlignment','bottom');
    end
    if ~isnan(V_cr(i))
        xline(V_cr(i), '--', 'Color', colors{i}, ...
              'Label', sprintf('V_{cr} = %.1f m/s', V_cr(i)), ...
              'LabelVerticalAlignment','bottom');
    end
end
xlabel('V  [m/s]');
ylabel('r / \delta   [rad/s per rad]');
title('Yaw rate gain');
legend('Location','northwest');

% =====================================================================
% PLOT 2 : SIDESLIP GAIN  beta/delta
% =====================================================================
figure('Color','w','Name','Sideslip gain');
hold on; grid on; box on;
for i = 1:3
    plot(V_vec, rad2deg(G_beta(i,:)), 'LineWidth', 1.8, 'Color', colors{i}, ...
         'DisplayName', labels{i});
end
for i = 1:3
    if ~isnan(V_ch(i))
        xline(V_ch(i), '--', 'Color', colors{i});
    end
    if ~isnan(V_cr(i))
        xline(V_cr(i), '--', 'Color', colors{i});
    end
end
xlabel('V  [m/s]');
ylabel('\beta / \delta   [deg per rad]');
title('Sideslip gain');
legend('Location','northwest');

% =====================================================================
% PLOT 3 : CURVATURE RESPONSE  (1/R)/delta
% =====================================================================
figure('Color','w','Name','Curvature response');
hold on; grid on; box on;
for i = 1:3
    plot(V_vec, G_curv(i,:), 'LineWidth', 1.8, 'Color', colors{i}, ...
         'DisplayName', labels{i});
end
for i = 1:3
    if ~isnan(V_ch(i))
        xline(V_ch(i), '--', 'Color', colors{i});
    end
    if ~isnan(V_cr(i))
        xline(V_cr(i), '--', 'Color', colors{i});
    end
end
xlabel('V  [m/s]');
ylabel('(1/R) / \delta   [1/m per rad]');
title('Curvature response');
legend('Location','northwest');

% =====================================================================
% PLOT 4 : LATERAL ACCELERATION RESPONSE  (V^2/R)/delta
% =====================================================================
figure('Color','w','Name','Lateral acceleration response');
hold on; grid on; box on;
for i = 1:3
    plot(V_vec, G_ay(i,:), 'LineWidth', 1.8, 'Color', colors{i}, ...
         'DisplayName', labels{i});
end
for i = 1:3
    if ~isnan(V_ch(i))
        xline(V_ch(i), '--', 'Color', colors{i});
    end
    if ~isnan(V_cr(i))
        xline(V_cr(i), '--', 'Color', colors{i});
    end
end
xlabel('V  [m/s]');
ylabel('(V^2/R) / \delta   [m/s^2 per rad]');
title('Lateral acceleration response');
legend('Location','northwest');

% =====================================================================
% SUMMARY
% =====================================================================
fprintf('\n  Veh    K [rad/(m/s^2)]      V_ch or V_cr\n');
fprintf('  ---    ----------------     ---------------------\n');
for i = 1:3
    if ~isnan(V_ch(i))
        fprintf('  %s     %+0.6f          V_ch = %6.2f m/s\n', ...
                labels{i}, K(i), V_ch(i));
    elseif ~isnan(V_cr(i))
        fprintf('  %s     %+0.6f          V_cr = %6.2f m/s\n', ...
                labels{i}, K(i), V_cr(i));
    else
        fprintf('  %s     %+0.6f          neutral (no finite V)\n', ...
                labels{i}, K(i));
    end
end
fprintf('\n');