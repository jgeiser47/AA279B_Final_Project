% Joshua Geiser
clear all; close all; clc;

% MJD_0 = cal_to_MJD(2030, 1, 1, 0, 0, 0)
% Earth_rv_0 = get_planet_rv(3, MJD_0)
% Mars_rv_0 = get_planet_rv(4, MJD_0)

% Get contour lines data by running Lambert's algorithm across search space
FLAG_use_stored_data = 1;
if FLAG_use_stored_data
    load('lambert_data.mat');
else
    s = get_contour_data();
end
plot_contour_blank(s)

% Get info on optimal launch/arrival date and plot a marker on contour plot
[min_val, min_index] = min(s.Z(:));
[min_row, min_col] = ind2sub(size(s.Z), min_index);
min_dep_MJD = s.X(min_row,min_col);
min_arr_MJD = s.Y(min_row,min_col);
min_dep_cal = datetime(min_dep_MJD, 'convertfrom', 'modifiedjuliandate');
min_arr_cal = datetime(min_arr_MJD, 'convertfrom', 'modifiedjuliandate');
plot(gca, toDateNum(min_dep_MJD), toDateNum(min_arr_MJD), 'r*');

% Get launch period info and plot as a line on contour plot
logical_1 = (s.Y == min_arr_MJD);
logical_2 = (s.Z <= 15.0);
logical_3 = (s.X <= 63725);
total_logical = and(and(logical_1, logical_2), logical_3);
plot(gca, toDateNum(s.X(total_logical)), toDateNum(s.Y(total_logical)), 'r-');

% Get more intermediate data on launch period
launch_period_MJDs_launch = s.X(total_logical);
launch_period_MJDs_arrival = zeros(size(launch_period_MJDs_launch)) + min_arr_MJD;
launch_window_open = datetime(launch_period_MJDs_launch(1), 'convertfrom', 'modifiedjuliandate');
launch_window_close = datetime(launch_period_MJDs_launch(end), 'convertfrom', 'modifiedjuliandate');
launch_period_dt = launch_period_MJDs_launch(end) - launch_period_MJDs_launch(1);

% Print relevant launch period data
fprintf('Launch Window Open:\t %s \n', launch_window_open);
fprintf('Launch Window Close:\t %s \n', launch_window_close);
fprintf('Launch Period Duration:\t %d [days] \n', launch_period_dt); 
fprintf('Mars Arrival Date:\t %s \n', min_arr_cal);

% RLA Contour Plot
figure(); hold on; grid on; axis equal;
levels = [0:10:360];
contour(toDateNum(s.X), toDateNum(s.Y), wrapTo360(rad2deg(s.RLA)), levels);
dateformat = 2;
datetick('x', dateformat);
datetick('y', dateformat);
xtickangle(45);
xlim([toDateNum(s.x_beg), toDateNum(s.x_end)]);
ylim([toDateNum(s.y_beg), toDateNum(s.y_end)]);
colorbar();
xlabel('Earth Departure Date'); ylabel('Mars Arrival Date');
title('RLA Contours for 2033 Earth-to-Mars launch');

% DLA Contour Plot
figure(); hold on; grid on; axis equal;
levels = [-90:5:90];
contour(toDateNum(s.X), toDateNum(s.Y), wrapTo180(rad2deg(s.DLA)), levels);
dateformat = 2;
datetick('x', dateformat);
datetick('y', dateformat);
xtickangle(45);
xlim([toDateNum(s.x_beg), toDateNum(s.x_end)]);
ylim([toDateNum(s.y_beg), toDateNum(s.y_end)]);
colorbar();
xlabel('Earth Departure Date'); ylabel('Mars Arrival Date');
title('DLA Contours for 2033 Earth-to-Mars launch');

% Launch Period RLA/DLA Plot
figure(); hold on; grid on; 
plot(toDateNum(s.X(total_logical)), wrapTo360(rad2deg(s.RLA(total_logical))), 'r-');
plot(toDateNum(s.X(total_logical)), wrapTo180(rad2deg(s.DLA(total_logical))), 'b-');
dateformat = 2;
datetick('x', dateformat);
xtickangle(45);
xlim([toDateNum(launch_period_MJDs_launch(1)), toDateNum(launch_period_MJDs_launch(end))]);
legend('RLA', 'DLA');
xlabel('Earth Departure Date'); ylabel('RLA/DLA [deg]');
title('RLA and DLA throughout Launch Period');
