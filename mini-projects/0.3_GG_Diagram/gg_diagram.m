function [ax_bnd, ay_bnd] = ggv_diagram(V, mu, aero_enabled)

% Vehicle parameters from linear_bicycle_model()
m  = 1500;       % vehicle mass [kg]
Iz = 2500;       % yaw moment of inertia [kg m^2]
a  = 1.2;        % CG -> front axle [m]
b  = 1.4;        % CG -> rear axle [m]
L  = a + b;      % wheelbase [m]

% Additional vehicle parameters
g  = 9.81;       % gravity [m/s^2]
h  = 0.5;        % CG height [m]
tf = 1.6;        % front track width [m]
tr = 1.6;        % rear track width [m]

% Aero
rho = 1.225;
ClA = 1.0;
aero_distribution_front = 0.5;

%% Aero load

if aero_enabled
    F_aero = 0.5*rho*ClA*V^2;
    F_aero_f = aero_distribution_front*F_aero;
    F_aero_r = (1-aero_distribution_front)*F_aero;
else
    F_aero_f = 0;
    F_aero_r = 0;
end

%% Acceleration direction

theta = linspace(0,2*pi,361);
R_max = zeros(size(theta));

%% Search maximum acceleration magnitude

for j = 1:length(theta)

    R_low  = 0;
    R_high = 30;

    for k = 1:60

        R = 0.5*(R_low + R_high);

        ax = R*cos(theta(j));
        ay = R*sin(theta(j));

        %% Total vehicle forces

        Fx = m*ax;
        Fy = m*ay;

        %% Longitudinal load transfer

        Fzf = (m*g*b - m*h*ax)/L + F_aero_f;
        Fzr = (m*g*a + m*h*ax)/L + F_aero_r;

        %% Lateral load transfer

        dFzf = m*abs(ay)*h*b/(L*tf);
        dFzr = m*abs(ay)*h*a/(L*tr);

        %% Individual wheel normal loads

        Fz_fl = Fzf/2 - dFzf;
        Fz_fr = Fzf/2 + dFzf;

        Fz_rl = Fzr/2 - dFzr;
        Fz_rr = Fzr/2 + dFzr;

        %% Bicycle-model lateral force distribution

        Fyf = (b/L)*Fy;
        Fyr = (a/L)*Fy;

        %% Rear-wheel drive

        Fxf = 0;
        Fxr = Fx;

        %% Equal left/right force distribution

        Fyf_l = Fyf/2;
        Fyf_r = Fyf/2;

        Fyr_l = Fyr/2;
        Fyr_r = Fyr/2;

        Fxr_l = Fxr/2;
        Fxr_r = Fxr/2;

        %% Friction ellipses

        front_left = (Fyf_l/(mu*Fz_fl))^2;

        front_right = (Fyf_r/(mu*Fz_fr))^2;

        rear_left = (Fxr_l/(mu*Fz_rl))^2 + (Fyr_l/(mu*Fz_rl))^2;

        rear_right = (Fxr_r/(mu*Fz_rr))^2 + (Fyr_r/(mu*Fz_rr))^2;

        %% Feasibility

        feasible = Fz_fl > 0 && Fz_fr > 0 && Fz_rl > 0 && ...
            Fz_rr > 0 && front_left  <= 1 && front_right <= 1 && ...
            rear_left   <= 1 && rear_right  <= 1;

        %% Binary search

        if feasible
            R_low = R;
        else
            R_high = R;
        end

    end

    R_max(j) = R_low;

end

%% GG boundary

ax_bnd = R_max .* cos(theta);
ay_bnd = R_max .* sin(theta);

plot(ax_bnd, ay_bnd)
axis equal
grid on
xlabel('a_x [m/s^2]')
ylabel('a_y [m/s^2]')

end