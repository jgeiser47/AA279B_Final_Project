% Purpose: Compare ephemeris and circular models of Mars' orbit to try to
% find a circular orbit that is reasonably close to the actual orbit. 
clear all; close all; clc;

% Assumes heliocentric orbit
mu_Sun = 1.3271244004193938e11;

% Mars at arrival epoch
MJD_Mars_arr  = cal_to_MJD(2034, 1, 28, 0, 0, 0);
rv_Mars_old  = get_planet_rv(4, MJD_Mars_arr);

a = 227953237.299835; 
e = 0.0934209570652977;
i = 0.0322348483513284;
RAAN = 0.863237189337131;
AOP = -1.27848943960252;
M = 1.06623161056081;

a = 227953016; n = sqrt(mu_Sun/(a^3));
e = 0;
M0 = M; T0 = MJD_Mars_arr;
rv_Mars_new = oe_to_rv(a, e, i, RAAN, AOP, M);

figure(); hold on; grid on; axis equal;
plot3(0,0,0,'y*');
plot3(rv_Mars_old(1), rv_Mars_old(2), rv_Mars_old(3), 'g*');
plot3(rv_Mars_new(1), rv_Mars_new(2), rv_Mars_new(3), 'k*');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
legend({'Sun', 'Ephemeris', 'Circular'}, 'Location', 'Northwest');
title('Difference between Mars Ephemeris/Circular Models on 01/28/2034');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MJD_0 = cal_to_MJD(2033, 3, 25, 0, 0, 0);
MJD_f = cal_to_MJD(2035, 1, 28, 0, 0, 0);
MJD_arr = [MJD_0 : 1 : MJD_f];

rv_Mars_old_vec = zeros(length(MJD_arr), 6);
rv_Mars_new_vec = zeros(length(MJD_arr), 6);

for ii = 1:length(MJD_arr)
    MJD_i = MJD_arr(ii);
    rv_Mars_old_vec(ii,:) = get_planet_rv(4, MJD_i)';
    
    Mi = M0 + n * ((MJD_i - T0)*86400);
    rv_Mars_new_vec(ii,:) = oe_to_rv(a, e, i, RAAN, AOP, Mi)';
    
end

figure(); hold on; grid on; axis equal;
plot3(0,0,0, 'y*');
plot3(rv_Mars_old_vec(:,1), rv_Mars_old_vec(:,2), rv_Mars_old_vec(:,3), 'g');
plot3(rv_Mars_new_vec(:,1), rv_Mars_new_vec(:,2), rv_Mars_new_vec(:,3), 'k');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
legend('Sun', 'Ephemeris', 'Circular');
title('Difference between Mars Ephemeris/Circular Models over ~1 Orbit');