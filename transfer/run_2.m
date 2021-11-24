% Joshua Geiser
clear all; close all; clc;

CONST = struct();

% Gravitational parameters
CONST.mu_Sun = 1.3271244004193938e11;
CONST.mu_Earth = 3.986004418e5;
CONST.mu_Mars = 42828.3719;
CONST.MJD_0 = cal_to_MJD(2033, 3, 25, 0, 0, 0);

% Earth at depature epoch
MJD_Earth_dep = cal_to_MJD(2033, 3, 25, 0, 0, 0);
rv_Earth_dep = get_planet_rv(3, MJD_Earth_dep);
r_Earth_dep = rv_Earth_dep(1:3);
v_Earth_dep = rv_Earth_dep(4:6);
CONST.r_Earth_dep = r_Earth_dep;
CONST.v_Earth_dep = v_Earth_dep;

% Mars at arrival epoch
MJD_Mars_arr  = cal_to_MJD(2034, 1, 28, 0, 0, 0);
rv_Mars_arr  = get_planet_rv(4, MJD_Mars_arr);
r_Mars_arr = rv_Mars_arr(1:3);
v_Mars_arr = rv_Mars_arr(4:6);
CONST.r_Mars_arr = r_Mars_arr;
CONST.v_Mars_arr = v_Mars_arr;

% Time of flight of trajectory
TOF = (MJD_Mars_arr - MJD_Earth_dep) * 86400;
CONST.TOF = TOF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From lambert solver solution
v0_guess = [4.74960464561663; -32.2010671899605; 0.0608576993798047]; 
plot_single_traj(CONST, v0_guess);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shooting method solution
x0 = v0_guess;
dt=TOF;
y_des = [r_Mars_arr(1)+10000; r_Mars_arr(2); r_Mars_arr(3)];
myfun = @run_traj_wrapper;
v0_converged = shooting(myfun,y_des,x0,dt,CONST);
plot_single_traj(CONST, v0_converged);

function plot_single_traj(CONST, v0)

r_Earth_dep = CONST.r_Earth_dep; TOF = CONST.TOF; MJD_Earth_dep = CONST.MJD_0;

x0 = [r_Earth_dep(1); r_Earth_dep(2); r_Earth_dep(3); v0(1); v0(2); v0(3)];
t_sim = [0 : 86400 : TOF];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_out, x_out] = ode113(@(t,x) calc_xdot_dim(t,x,CONST), t_sim, x0, options);

r_Earth_vec = zeros(length(t_sim),6);
r_Mars_vec = zeros(length(t_sim),6);
for i = 1:length(t_sim)
    MJD_i = MJD_Earth_dep + t_sim(i)/86400;
    r_Earth_vec(i,:) = get_planet_rv(3, MJD_i)';
    r_Mars_vec(i,:) = get_planet_rv(4, MJD_i)';
end

figure(); hold on; grid on; axis equal;
plot3(r_Earth_vec(:,1), r_Earth_vec(:,2), r_Earth_vec(:,3), 'b:', 'Linewidth', 2);
plot3(r_Mars_vec(:,1), r_Mars_vec(:,2), r_Mars_vec(:,3), 'r:', 'Linewidth', 2);
plot3(x_out(:,1), x_out(:,2), x_out(:,3), 'k-', 'Linewidth', 2);

error = norm(r_Mars_vec(end,1:3) - x_out(end,1:3))
end

function x0_out = shooting(f,y_des,x0,TOF,CONST)

% Shooting method parameters
conv_tol = 1e-6;
max_iters = 100;
verbose = 1;

iter = 0;
while (iter < max_iters)
    
    % Increment iteration and calculate current error
    iter = iter + 1;
    y_out = f(x0, TOF, CONST);
    delta_y = y_out - y_des;
    
    % Print iteration/error data
    if verbose
        disp(['Iteration: ', num2str(iter), sprintf('\t'), ...
              'Error Norm: ', num2str(norm(delta_y))]);
    end
    
    % End if converged
    if (norm(delta_y) < conv_tol)
        break
    end
    
    % If not converged, calculate Jacobian and adjust input vector
    J = calc_J(f,x0,TOF,CONST);
    x0 = x0 - pinv(J) * (delta_y);
end
x0_out = x0;
end

function J = calc_J(f,x0,dt,CONST)

N_in = length(x0);          % length of input vector
N_out = 3;                  % length of output vector
J = zeros(N_in, N_out);     % Jacobian matrix

% Loop through each variable in input vector
for i = 1:N_in
    
    % Base the perturbation size on the magnitude of the input variable
    perturbation = zeros(N_in,1);
    x_var_curr = x0(i);
    order_of_x_var_curr = floor(log(abs(x_var_curr))./log(10));
    delta = 10^(order_of_x_var_curr-4);
    perturbation(i) = delta;
    
    % Use finite differencing to calculate variations in output variables
    y_out_high = f(x0+perturbation, dt, CONST);
    y_out_low  = f(x0-perturbation, dt, CONST);
    J(:,i) = (y_out_high - y_out_low) ./ (2*delta);
end
end

function x_f = run_traj_wrapper(x0,TOF,CONST)

% Setup of initial conditions
x0_full = [CONST.r_Earth_dep; x0];

% Run sim
t_sim = [0 TOF];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_out, x_out] = ode113(@(t,x) calc_xdot_dim(t,x,CONST), t_sim, x0_full, options); 

% Only output we care about is final (X,Y) position
x_f = [x_out(end,1); x_out(end,2); x_out(end,3)];
end

function xdot = calc_xdot_dim(t,x,CONST)
% Author: Josh Geiser
% Inputs: t     - current timestep
%         x     - current state vector
% Output: xdot  - time derivate of state vector at current timestep

% Constants
mu_Sun = CONST.mu_Sun;
mu_Earth = CONST.mu_Earth;
mu_Mars = CONST.mu_Mars;
MJD_0 = CONST.MJD_0;

% Unpack state vector
X = x(1); Y = x(2); Z = x(3); VX = x(4); VY = x(5); VZ = x(6);

% Intermediate calculations
v = [VX; VY; VZ];
r = [X; Y; Z];
MJD = MJD_0 + t/86400;
rv_Earth = get_planet_rv(3, MJD);
r_Earth = rv_Earth(1:3);
rv_Mars = get_planet_rv(4, MJD);
r_Mars = rv_Mars(1:3);
r2 = r - r_Mars;

% State space form
xdot = zeros(6,1);
xdot(1:3) = v;
xdot(4:6) = ((-mu_Sun/(norm(r)^3)) .* r) - ((-mu_Mars/(norm(r2)^3)) .* r2);
end
