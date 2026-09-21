function [A,B] = linear_bicycle_model(V)
% BICYCLE_MODEL_MATRICES
% Linear 2-DOF bicycle model — state-space matrices A and B.
%
%   [A,B] = linear_bicycle_model(V)
%
%    arguments to this function :  p, V          →  parameters that define A and B
%    state vector               :  X = [v; r]    →  appears as X in dX = A·X + B·δ
%    input                      :  δ             →  appears as the input in B·δ
%
%
% =========================================================================
% DIMENSIONS OF A AND B
% =========================================================================
%   A11 : (N/rad) / (kg*m/s)            = 1/s        (rad is dimensionless)
%   A12 : (m*N/rad)/(kg*m/s) - m/s      = m/s
%   A21 : (m*N/rad) / (kg*m^2 * m/s)    = 1/(m*s)
%   A22 : (m^2*N/rad)/(kg*m^2 * m/s)    = 1/s
%   B1  : (N/rad)/kg                    = m/(s^2*rad)
%   B2  : (m*N/rad)/(kg*m^2)            = 1/(s^2*rad)
%
%  With X = [m/s ; rad/s] and delta = rad, both rows of dX carry the
%  correct units: m/s^2 and rad/s^2 respectively.
%
% =========================================================================
% PARAMETER STRUCT
% =========================================================================
%   m   : total mass                             [kg]
%   Iz  : yaw moment of inertia about CG          [kg*m^2]
%   a   : CG to front axle                        [m]
%   b   : CG to rear axle                         [m]
%   CF  : front axle cornering stiffness          [N/rad]
%   CR  : rear axle cornering stiffness           [N/rad]
%
% =========================================================================

m   = 1500;
Iz  = 2500;
a   = 1.2;
b   = 1.4;
CF  = 80000;
CR  = 90000;

% --- Stability derivatives --------------------------------------------
Y_beta  = CF + CR;
Y_r     = (a*CF - b*CR)/V;
Y_delta = -CF;
N_beta  = a*CF - b*CR;
N_r     = (a^2*CF + b^2*CR)/V;
N_delta = -a*CF;

% --- State-space matrices ---------------------------------------------
A = [ Y_beta/(m*V),    Y_r/m - V ;
      N_beta/(Iz*V),   N_r/Iz    ];

B = [ Y_delta/m ;
      N_delta/Iz ];

end