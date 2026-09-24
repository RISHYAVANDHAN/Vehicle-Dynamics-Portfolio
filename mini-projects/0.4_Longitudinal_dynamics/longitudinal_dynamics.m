% Longitudinal dynamics solver
% Point mass Vehicle for the longitudinal dynamics solver
% We are providing the point mass states x = [v:s], v = vehilce velocity, s = vehicle position 
% and learning how the vehicle's motion is in the longitudinal axis i.e
% forward backward.
% F_tractive/drive = F_rr + F_aero + F_accelration
% F= ma = F_accel = F_drive - (F_rr + F_aero)
% x = [v;s] the vehicle is moving at 30m/s and positioned at s = 3m 

% same bicycle model vehicle parameters
m = 1500;           % mass of the point mass vehicle 
Crr = 0.01;         % rolling resistance of a passenger car tire
Cd = 0.3;           % Aero drag
rho = 1.224;        % Air density
A = 2;              % Frontal area = 2m^2
F_drive = 3500;     % Tractive force from the engine

t_span = [0 60]; 
x0 = [0; 0]; 
v_target = 100 / 3.6; % m/s to km/hr for a 0-100 km/hr sprint time calculation

% MATLAB Hint: Solver configuration settings
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%% 1. Baseline 0–100 km/h Run & Plotting
Cd = 0.3;
[t_base, X_base] = ode45(@(t,x) longitudinal_dynamics(t, x, F_drive, m, Crr, Cd, rho, A), t_span, x0, opts);
v_base = X_base(:, 1);

% Recompute acceleration array for plotting deliverables
a_base = zeros(length(t_base), 1);
for i = 1:length(t_base)
    xdot = longitudinal_dynamics(t_base(i), X_base(i,:)', F_drive, m, Crr, Cd, rho, A);
    a_base(i) = xdot(1);
end

% MATLAB Hint: Interpolate exact sprint timeline
t_0_100 = interp1(v_base, t_base, v_target);

% Plot execution
figure; 
yyaxis left; 
plot(t_base, v_base * 3.6, 'b-', 'LineWidth', 2); 
ylabel('Velocity (km/h)');
yyaxis right; 
plot(t_base, a_base, 'r-', 'LineWidth', 2); 
ylabel('Acceleration (m/s^2)');
xlabel('Time (s)'); 
title(['0-100 km/h Sprint Time: ', num2str(t_0_100), ' s']); 
grid on;

%% 2. Zero-Drag Analytical Validation
Cd_zero = 0;
[t_zero, X_zero] = ode45(@(t,x) longitudinal_dynamics(t, x, F_drive, m, Crr, Cd_zero, rho, A), t_span, x0, opts);

a_const = (F_drive - (Crr * m * 9.81)) / m;
v_analyt = a_const .* t_zero;
s_analyt = 0.5 .* a_const .* (t_zero.^2);

v_err = max(abs(X_zero(:,1) - v_analyt));
s_err = max(abs(X_zero(:,2) - s_analyt));
fprintf('Validation -> Velocity Error: %.4e m/s | Displacement Error: %.4e m\n\n', v_err, s_err);

%% 3. Sensitivity Analysis Matrix
Cds = [0.25, 0.30, 0.35];
fprintf('%-8s | %-15s | %-15s\n', 'Cd', '0-100 Time (s)', 'V_max (km/h)');
for Cd_val = Cds
    [t_s, X_s] = ode45(@(t,x) longitudinal_dynamics(t, x, F_drive, m, Crr, Cd_val, rho, A), t_span, x0, opts);
    t_sprint = interp1(X_s(:,1), t_s, v_target);
    v_max = sqrt((2 * (F_drive - (Crr * m * 9.81))) / (rho * Cd_val * A)) * 3.6;
    fprintf('%.2f     | %.2f           | %.1f\n', Cd_val, t_sprint, v_max);
end

%% --- ODE Core Function Definition ---
function xdot = longitudinal_dynamics(t, x, F_drive, m, Crr, Cd, rho, A)
    v = x(1);
    F_drag = 0.5 * rho * Cd * A * (v^2);
    F_roll = Crr * m * 9.81;
    a = (F_drive - F_drag - F_roll) / m;
    xdot = [a; v];
end
