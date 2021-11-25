% Run test NRHO plotting
clear all; close all; clc;

load('nominal_north_NRHO.mat');

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

% Dimensionalized quantities
a_Mars = 227953016; n_Mars = sqrt(CONST.mu_Sun/(a_Mars^3));

X_nrho = X_nrho(:,1:3);
X_nrho = X_nrho .* a_Mars;

figure(); hold on; grid on; axis equal;
plot3(0,0,0,'y*');
plot3(a_Mars,0,0,'r*');
plot3(X_nrho(:,1), X_nrho(:,2), X_nrho(:,3), 'k-');