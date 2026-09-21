% LEARN_BICYCLE_MODEL
% Guided derivation of the linear 2-DOF bicycle model.
%
% Run this file top to bottom. It will:
%   1. Build the equations of motion symbolically.
%   2. Print each intermediate expression (slip angles, tyre forces,
%      total force, yaw moment, state derivatives).
%   3. Extract A and B by differentiation.
%   4. Substitute a concrete vehicle at 25 m/s and show the numeric
%      result with the units of each element.
%   5. Cross-check against bicycle_model_matrices.m.
%
% Sign convention:
%   alpha_F = (v + a*r)/V - delta
%   alpha_R = (v - b*r)/V
%   Y_F = CF*alpha_F,  Y_R = CR*alpha_R
%   Y   = Y_F + Y_R,   N = a*Y_F - b*Y_R
%   m*(v_dot + V*r) = Y,   Iz*r_dot = N
%
% The printed expressions should match, line by line, what you wrote by hand.

clear; clc;
fprintf('\n====================================================\n');
fprintf('  LINEAR 2-DOF BICYCLE MODEL — GUIDED DERIVATION\n');
fprintf('====================================================\n\n');

%% --- 0. Parameters for the concrete example ---------------------------
% Define these up front so they can be substituted at the end.
p.m  = 1500;      % kg
p.Iz = 2500;      % kg*m^2
p.a  = 1.2;       % m
p.b  = 1.4;       % m
p.CF = 80000;     % N/rad
p.CR = 90000;     % N/rad
V_val = 25;       % m/s

fprintf('Concrete vehicle for the numeric example:\n');
fprintf('  m   = %.0f kg\n',  p.m);
fprintf('  Iz  = %.0f kg*m^2\n', p.Iz);
fprintf('  a   = %.2f m\n',  p.a);
fprintf('  b   = %.2f m\n',  p.b);
fprintf('  CF  = %.0f N/rad\n', p.CF);
fprintf('  CR  = %.0f N/rad\n', p.CR);
fprintf('  V   = %.0f m/s\n\n', V_val);

%% --- 1. Symbolic variables -------------------------------------------
syms m Iz a b CF CR V v r delta v_dot r_dot real

fprintf('----------------------------------------------------\n');
fprintf(' STEP 1.  Slip angles  \n');
fprintf('----------------------------------------------------\n');
alpha_F = (v + a*r)/V - delta;
alpha_R = (v - b*r)/V;
fprintf('  alpha_F = '); disp(alpha_F);
fprintf('  alpha_R = '); disp(alpha_R);

fprintf('----------------------------------------------------\n');
fprintf(' STEP 2.  Linear tyre forces\n');
fprintf('----------------------------------------------------\n');
Y_F = CF * alpha_F;
Y_R = CR * alpha_R;
Y_F_expanded = expand(Y_F);
Y_R_expanded = expand(Y_R);
fprintf('  Y_F = CF*alpha_F = '); disp(Y_F_expanded);
fprintf('  Y_R = CR*alpha_R = '); disp(Y_R_expanded);

fprintf('----------------------------------------------------\n');
fprintf(' STEP 3.  Total lateral force  \n');
fprintf('----------------------------------------------------\n');
Y = Y_F + Y_R;
Y = collect(Y, [v, r, delta]);
fprintf('  Y = Y_F + Y_R = '); disp(Y);

fprintf('----------------------------------------------------\n');
fprintf(' STEP 4.  Total yaw moment about CG  \n');
fprintf('----------------------------------------------------\n');
N = a*Y_F - b*Y_R;
N = collect(N, [v, r, delta]);
fprintf('  N = a*Y_F - b*Y_R = '); disp(N);

fprintf('----------------------------------------------------\n');
fprintf(' STEP 5.  Equations of motion  \n');
fprintf('----------------------------------------------------\n');
eq_lat = m*(v_dot + V*r) == Y;
eq_yaw = Iz*r_dot        == N;
fprintf('  Lateral : '); disp(eq_lat);
fprintf('  Yaw     : '); disp(eq_yaw);

fprintf('----------------------------------------------------\n');
fprintf(' STEP 6.  Solve for v_dot and r_dot\n');
fprintf('----------------------------------------------------\n');
S = solve([eq_lat, eq_yaw], [v_dot, r_dot]);
v_dot_expr = simplify(S.v_dot);
r_dot_expr = simplify(S.r_dot);
fprintf('  v_dot = '); disp(v_dot_expr);
fprintf('  r_dot = '); disp(r_dot_expr);

fprintf('----------------------------------------------------\n');
fprintf(' STEP 7.  Read off A and B\n');
fprintf('         dX = A*X + B*delta,  X = [v; r]\n');
fprintf('----------------------------------------------------\n');
A_sym = simplify(jacobian([v_dot_expr; r_dot_expr], [v, r]));
B_sym = simplify(jacobian([v_dot_expr; r_dot_expr], delta));
fprintf('  Symbolic A =\n'); disp(A_sym);
fprintf('  Symbolic B =\n'); disp(B_sym);

fprintf('----------------------------------------------------\n');
fprintf(' STEP 8.  Stability derivatives \n');
fprintf('----------------------------------------------------\n');
Y_beta  = simplify(diff(Y, v) * V);       % dY/dbeta, beta = v/V
Y_r_sym = simplify(diff(Y, r));
Y_d_sym = simplify(diff(Y, delta));
N_beta  = simplify(diff(N, v) * V);
N_r_sym = simplify(diff(N, r));
N_d_sym = simplify(diff(N, delta));
fprintf('  Y_beta  = CF + CR                 -> '); disp(Y_beta);
fprintf('  Y_r     = (a*CF - b*CR)/V         -> '); disp(Y_r_sym);
fprintf('  Y_delta = -CF                     -> '); disp(Y_d_sym);
fprintf('  N_beta  = a*CF - b*CR             -> '); disp(N_beta);
fprintf('  N_r     = (a^2*CF + b^2*CR)/V     -> '); disp(N_r_sym);
fprintf('  N_delta = -a*CF                   -> '); disp(N_d_sym);

%% --- 2. Numeric substitution -----------------------------------------
fprintf('\n----------------------------------------------------\n');
fprintf(' STEP 9.  Numeric A and B for the concrete vehicle\n');
fprintf('         V = %.1f m/s\n', V_val);
fprintf('----------------------------------------------------\n');

subs_vals  = {m, Iz, a, b, CF, CR, V};
subs_nums  = {p.m, p.Iz, p.a, p.b, p.CF, p.CR, V_val};

A_num = double(subs(A_sym, subs_vals, subs_nums));
B_num = double(subs(B_sym, subs_vals, subs_nums));

fprintf('  A =\n');
fprintf('    [ %+10.4f   %+10.4f ]\n', A_num(1,1), A_num(1,2));
fprintf('    [ %+10.4f   %+10.4f ]\n', A_num(2,1), A_num(2,2));
fprintf('  B =\n');
fprintf('    [ %+10.4f ]\n', B_num(1));
fprintf('    [ %+10.4f ]\n', B_num(2));

fprintf('\n  Element-by-element with units:\n');
fprintf('    A11 = %+8.4f   1/s\n',   A_num(1,1));
fprintf('    A12 = %+8.4f   m/s\n',   A_num(1,2));
fprintf('    A21 = %+8.4f   1/(m*s)\n', A_num(2,1));
fprintf('    A22 = %+8.4f   1/s\n',   A_num(2,2));
fprintf('    B1  = %+8.4f   m/(s^2*rad)\n', B_num(1));
fprintf('    B2  = %+8.4f   1/(s^2*rad)\n', B_num(2));

%% --- 3. Physical sanity checks ---------------------------------------
fprintf('\n----------------------------------------------------\n');
fprintf(' STEP 10. Physical sanity checks\n');
fprintf('----------------------------------------------------\n');

% Understeer gradient K = m/L^2 * (a/CR - b/CF)
L = p.a + p.b;
K = p.m / L^2 * (p.a/p.CR - p.b/p.CF);
fprintf('  Understeer gradient K = %+0.5f rad/(m/s^2)\n', K);
if K > 0
    fprintf('    -> K > 0 : understeer\n');
    fprintf('       V_char = %.1f m/s\n', sqrt(1/K));
elseif K < 0
    fprintf('    -> K < 0 : oversteer\n');
    fprintf('       V_crit = %.1f m/s\n', sqrt(-1/K));
else
    fprintf('    -> K = 0 : neutral steer\n');
end

% Eigenvalues of A -> stability and mode character
lambda = eig(A_num);
fprintf('\n  Eigenvalues of A at V = %.1f m/s:\n', V_val);
for k = 1:numel(lambda)
    fprintf('    lambda_%d = %+0.4f %+0.4fi\n', ...
            k, real(lambda(k)), imag(lambda(k)));
end
if all(real(lambda) < 0)
    fprintf('    -> all Re(lambda) < 0 : stable\n');
else
    fprintf('    -> at least one Re(lambda) > 0 : unstable\n');
end

%% --- 4. Cross-check against the production function ------------------
fprintf('\n----------------------------------------------------\n');
fprintf(' STEP 11. Cross-check with bicycle_model_matrices.m\n');
fprintf('----------------------------------------------------\n');
[A_num2, B_num2] = bicycle_model_matrices(V_val);
errA = max(abs(A_num(:) - A_num2(:)));
errB = max(abs(B_num(:) - B_num2(:)));
fprintf('  max |A_sym - A_num| = %.3e\n', errA);
fprintf('  max |B_sym - B_num| = %.3e\n', errB);
assert(errA < 1e-10 && errB < 1e-10, ...
       'Mismatch between symbolic and numeric implementations.');
fprintf('  -> matrices match. Derivation is consistent.\n\n');