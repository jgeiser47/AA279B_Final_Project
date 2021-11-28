% Joshua Geiser
% Launch Window Calculations 
clear all; close all; clc;

% Load lambert data structure 's'
load('lambert_data.mat');

% Nominal launch and landing MJDs
launch_MJD = cal_to_MJD(2033, 3, 28, 0, 0, 0);
land_MJD = cal_to_MJD(2034, 2, 15, 0, 0, 0);

% Range safety limits for KSC in degrees
min_lim = 70;
max_lim = 115;

% Get azimuth value at each launch epoch
[MJD_arr, beta_arr] = run_launch_date(launch_MJD, land_MJD, s);

% Find when beta is within range safety limits
log_inds = and((beta_arr>min_lim), (beta_arr<max_lim));

% Make a plot
figure(); hold on; grid on;
patch([-100 100 100 -100], [min_lim min_lim -180 -180], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
patch([-100 100 100 -100], [max_lim max_lim 180 180], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
plot((MJD_arr-launch_MJD)*24, beta_arr, 'ro', 'Markersize', 2);
plot((MJD_arr(log_inds)-launch_MJD)*24, beta_arr(log_inds), 'bo', 'Markersize', 2);
xlim([0 25]), ylim([0 180]);
xlabel('UTC Time Past Midnight [hours]');
ylabel('Launch Azimuth [deg]');
title('Launch Azimuth vs Launch Time for 03/28/2033 Launch Date');

valid_betas = MJD_arr(log_inds);
temp_valid_betas = abs(diff(valid_betas));
[~,split_ind] = max(temp_valid_betas);
LW_1 = valid_betas(1:split_ind);
LW_2 = valid_betas(split_ind+1:end);

LW_1_open_UTC = datetime(LW_1(1), 'convertfrom', 'modifiedjuliandate');
LW_1_clos_UTC = datetime(LW_1(end), 'convertfrom', 'modifiedjuliandate');
LW_1_open_EDT = datetime(LW_1(1)-(4/24), 'convertfrom', 'modifiedjuliandate');
LW_1_clos_EDT = datetime(LW_1(end)-(4/24), 'convertfrom', 'modifiedjuliandate');
LW_2_open_UTC = datetime(LW_2(1), 'convertfrom', 'modifiedjuliandate');
LW_2_clos_UTC = datetime(LW_2(end), 'convertfrom', 'modifiedjuliandate');
LW_2_open_EDT = datetime(LW_2(1)-(4/24), 'convertfrom', 'modifiedjuliandate');
LW_2_clos_EDT = datetime(LW_2(end)-(4/24), 'convertfrom', 'modifiedjuliandate');

fprintf('Launch Window 1 Open UTC: \t %s \n', LW_1_open_UTC);
fprintf('Launch Window 1 Close UTC:\t %s \n', LW_1_clos_UTC);
fprintf('Launch Window 1 Open EDT: \t %s \n', LW_1_open_EDT);
fprintf('Launch Window 1 Close EDT:\t %s \n', LW_1_clos_EDT);
fprintf('\n');
fprintf('Launch Window 2 Open UTC: \t %s \n', LW_2_open_UTC);
fprintf('Launch Window 2 Close UTC:\t %s \n', LW_2_clos_UTC);
fprintf('Launch Window 2 Open EDT: \t %s \n', LW_2_open_EDT);
fprintf('Launch Window 2 Close EDT:\t %s \n', LW_2_clos_EDT);

function [MJD_arr, beta_arr]  = run_launch_date(launch_date, land_MJD, s)

% First need to get RLA/DLA values (assumed constant over 1 launch date)
logical_1 = (s.X == launch_date);
logical_2 = (s.Y == land_MJD);
total_logical = and(logical_1, logical_2);
RLA = s.RLA(total_logical);
DLA = s.DLA(total_logical);

MJD_arr = (launch_date + ([3600*2.25 : 60 : 3600*26.175] ./86400))';
beta_arr = zeros(length(MJD_arr), 1);

for i = 1:length(MJD_arr)
    beta_arr(i) = get_beta(MJD_arr(i), RLA, DLA);
end
end


function beta = get_beta(launch_time, RLA, DLA)

% Reference values
MJD_ref = cal_to_MJD(2018, 01, 01, 00, 00, 00);
GMST_ref = 6.7066152 * (2*pi/24);
launch_lat = deg2rad(28.573469);
launch_lon = deg2rad(-80.651070); 
OM_EARTH = 0.0000729211585530; 

% Launch azimuth at launch time
alpha_l = wrapTo2Pi(launch_lon + GMST_ref + OM_EARTH * ((launch_time-MJD_ref) * 86400));
beta = atan2( sin(RLA-alpha_l) , cos(launch_lat)*tan(DLA) - sin(launch_lat)*cos(RLA-alpha_l) );
beta = rad2deg(beta);
if beta < 0
    beta = beta + 180;
end
end