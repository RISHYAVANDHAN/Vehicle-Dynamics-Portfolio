function xdot = longitudinal_dynamics(t,x,f_drive)
% Point mass Vehicle for the longitudinal dynamics solver
% We are providing the point mass states x = [v:s], v = vehilce velocity, s = vehicle position 
% and learning how the vehicle's motion is in the longitudinal axis i.e
% forward backward.
% F_tractive/drive = F_rr + F_aero + F_grade + F_accelration
% F= ma = F_accel = F_drive - (F_rr + F_aero + F_grade)
% same bicycle model vehicle parameters
m = 1500;           % mass of the point mass vehicle 
t = 60;             % runtime for the ode solver
v = 30;             % x = [v;s] the vehicle is moving at 30m/s and positioned at s = 3m 
F_drive = 50;       % driving force, say the powertrain produces _N
Crr = 0.01;         % rolling resistance of a passenger car tire
Cd = 0.3;           % Aero drag
rho = 1.224;        % Air density
A = 2;              % Frontal area = 2m^2

F_drag = 0.5 * rho * Cd * A * (v^2);
F_roll = Crr * m * 9.81;
a = (F_drive - F_drag - F_roll) / m;

xdot = [a;v];

end