% Joshua Geiser
clear all; close all; clc;

CONST = struct();

% Gravitational parameters
CONST.mu_Sun = 1.3271244004193938e11;
CONST.mu_Earth = 3.986004418e5;
CONST.mu_Mars = 42828.375816;
CONST.MJD_0 = cal_to_MJD(2033, 3, 25, 0, 0, 0);

% Earth at depature epoch
MJD_Earth_dep = cal_to_MJD(2033, 3, 25, 0, 0, 0);
rv_Earth_dep = get_Earth_rv(MJD_Earth_dep);
r_Earth_dep = rv_Earth_dep(1:3);
v_Earth_dep = rv_Earth_dep(4:6);
CONST.r_Earth_dep = r_Earth_dep;
CONST.v_Earth_dep = v_Earth_dep;

% Mars at arrival epoch
MJD_Mars_arr  = cal_to_MJD(2034, 1, 28, 0, 0, 0);
rv_Mars_arr  = get_Mars_rv(MJD_Mars_arr);
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
%plot_single_traj(CONST, v0_guess);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shooting method solution

a_Mars = 227953016; n_Mars = sqrt(CONST.mu_Sun/(a_Mars^3));
%y_des = [r_Mars_arr(1); r_Mars_arr(2); r_Mars_arr(3)+15000];
rv_des_syn_nondim = [1.000291664674990; 0; 0.005620030466350; 0; -0.000789269445422; 0];
rv_des_syn = rv_des_syn_nondim .* a_Mars; rv_des_syn(4:6) = rv_des_syn(4:6) .* n_Mars;
rv_des = syn_to_inert(MJD_Mars_arr, rv_des_syn')';

x0 = v0_guess;
dt=TOF;
myfun = @run_traj_wrapper;
v0_converged = shooting(myfun,rv_des(1:3),x0,dt,CONST);
plot_single_traj(CONST, v0_converged);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attempt at getting NRHO section working based on IC

% x0 = rv_des;
% t_sim = [TOF : 86400/10 : TOF+(86400*100)];
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% [t_NRHO, x_NRHO] = ode113(@(t,x) calc_xdot_dim(t,x,CONST), t_sim, x0, options);
% 
% MJDs_NRHO = (t_NRHO./86400) + CONST.MJD_0;
% x_syn_NRHO = inert_to_syn(MJDs_NRHO, x_NRHO(:,1:3));

% figure(); hold on; grid on; axis equal;
% plot3(gca,0,0,0,'y*');
% plot3(gca,a_Mars,0,0,'r*');
% plot3(gca,x_syn_NRHO(:,1), x_syn_NRHO(:,2), x_syn_NRHO(:,3), 'm-', 'Linewidth', 2);

function plot_single_traj(CONST, v0)

r_Earth_dep = CONST.r_Earth_dep; TOF = CONST.TOF; MJD_Earth_dep = CONST.MJD_0;

x0 = [r_Earth_dep(1); r_Earth_dep(2); r_Earth_dep(3); v0(1); v0(2); v0(3)];
t_sim = [0 : 86400 : TOF];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_out, x_out] = ode113(@(t,x) calc_xdot_dim(t,x,CONST), t_sim, x0, options);

r_Earth_vec = zeros(length(t_sim),6);
r_Mars_vec = zeros(length(t_sim),6);
MJDs_vec = zeros(length(t_sim),1);
for i = 1:length(t_sim)
    MJD_i = MJD_Earth_dep + t_sim(i)/86400;
    MJDs_vec(i) = MJD_i;
    r_Earth_vec(i,:) = get_Earth_rv(MJD_i)';
    r_Mars_vec(i,:) = get_Mars_rv(MJD_i)';
end

figure(); hold on; grid on; axis equal;
plot3(x_out(:,1), x_out(:,2), x_out(:,3), 'k-', 'Linewidth', 2);
plot3(r_Earth_vec(:,1), r_Earth_vec(:,2), r_Earth_vec(:,3), 'b:', 'Linewidth', 2);
plot3(r_Mars_vec(:,1), r_Mars_vec(:,2), r_Mars_vec(:,3), 'r:', 'Linewidth', 2);
plot3(0,0,0,'y*');
legend({'Spacecraft', 'Earth', 'Mars', 'Sun'}, 'Location', 'Northwest');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Transfer Trajectory in Heliocentric Intertial Frame');

%error = norm(r_Mars_vec(end,1:3) - x_out(end,1:3))

x_out_syn = inert_to_syn(MJDs_vec, x_out(:,1:3));
r_Earth_vec_syn = inert_to_syn(MJDs_vec, r_Earth_vec(:,1:3));
r_Mars_vec_syn = inert_to_syn(MJDs_vec, r_Mars_vec(:,1:3));

figure(); hold on; grid on; axis equal;
plot3(x_out_syn(:,1), x_out_syn(:,2), x_out_syn(:,3), 'k', 'Linewidth', 2);
plot3(r_Earth_vec_syn(:,1), r_Earth_vec_syn(:,2), r_Earth_vec_syn(:,3), 'b:', 'Linewidth', 2);
plot3(r_Mars_vec_syn(:,1), r_Mars_vec_syn(:,2), r_Mars_vec_syn(:,3), 'r*');
plot3(0,0,0,'y*');
legend({'Spacecraft', 'Earth', 'Mars', 'Sun'}, 'Location', 'Northeast');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Transfer Trajectory in Sun-Mars Rotating Frame');
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
rv_Earth = get_Earth_rv(MJD);
r_Earth = rv_Earth(1:3);
rv_Mars = get_Mars_rv(MJD);
r_Mars = rv_Mars(1:3);
r2 = r - r_Mars;
if norm(r2) < 9000
    disp('Warning: too close to Mars, may have trouble converging');
end

% State space form
xdot = zeros(6,1);
xdot(1:3) = v;
xdot(4:6) = ((-mu_Sun/(norm(r)^3)) .* r) - ((-mu_Mars/(norm(r2)^3)) .* r2);
end
