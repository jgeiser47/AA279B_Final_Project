% Joshua Geiser
clear all; close all; clc;

% Get contour lines data by running Lambert's algorithm across search space
FLAG_use_stored_data = 1;
if FLAG_use_stored_data
    load('lambert_data.mat');
else
    s = get_contour_data();
end
plot_contour_blank(s)
xlim([toDateNum(s.x_beg+32), toDateNum(s.x_end-32)]);
ylim([toDateNum(s.y_beg+32), toDateNum(s.y_end-32)]);

% Get info on optimal launch/arrival date and plot a marker on contour plot
[min_val, min_index] = min(s.Z(:));
[min_row, min_col] = ind2sub(size(s.Z), min_index);
min_dep_MJD = s.X(min_row,min_col);
min_arr_MJD = s.Y(min_row,min_col);
min_dep_cal = datetime(min_dep_MJD, 'convertfrom', 'modifiedjuliandate');
min_arr_cal = datetime(min_arr_MJD, 'convertfrom', 'modifiedjuliandate');

% Get launch period info and plot as a line on contour plot
logical_1 = (s.Y == min_arr_MJD);
logical_2 = (s.Z <= 15.0);
logical_3 = (s.X <= 63725);
total_logical = and(and(logical_1, logical_2), logical_3);

% Get more intermediate data on launch period
launch_period_MJDs_launch = s.X(total_logical);
launch_period_MJDs_arrival = zeros(size(launch_period_MJDs_launch)) + min_arr_MJD;
launch_window_open = datetime(launch_period_MJDs_launch(1), 'convertfrom', 'modifiedjuliandate');
launch_window_close = datetime(launch_period_MJDs_launch(end), 'convertfrom', 'modifiedjuliandate');
launch_period_dt = launch_period_MJDs_launch(end) - launch_period_MJDs_launch(1);

% If plotting specific launch period data
if 1
    xlim([toDateNum(s.x_beg+100), toDateNum(s.x_end-80)]);
    ylim([toDateNum(s.y_beg+180), toDateNum(s.y_end-100)]);
    plot(gca, toDateNum(s.X(total_logical)), toDateNum(s.Y(total_logical)), 'r-', 'Linewidth', 1);
    plot(gca, toDateNum(launch_period_MJDs_launch(1)), toDateNum(min_arr_MJD), 'rd', 'Markersize', 8);
    plot(gca, toDateNum(min_dep_MJD), toDateNum(min_arr_MJD), 'r*', 'Markersize', 8);
    plot(gca, toDateNum(launch_period_MJDs_launch(end)), toDateNum(min_arr_MJD), 'rx', 'Markersize', 8);
    title('Launch Period with C_3 Contours');
    legend({'Contours','Launch Period','Launch Period Open','Optimal Launch Date','Launch Period Close'}, 'Location', 'Northwest');
end

% Print relevant launch period data
fprintf('Launch Period Open:\t %s \n', launch_window_open);
fprintf('Launch Period Close:\t %s \n', launch_window_close);
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
xlim([toDateNum(s.x_beg+32), toDateNum(s.x_end-32)]);
ylim([toDateNum(s.y_beg+32), toDateNum(s.y_end-32)]);
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
xlim([toDateNum(s.x_beg+32), toDateNum(s.x_end-32)]);
ylim([toDateNum(s.y_beg+32), toDateNum(s.y_end-32)]);
colorbar();
xlabel('Earth Departure Date'); ylabel('Mars Arrival Date');
title('DLA Contours for 2033 Earth-to-Mars launch');

% Launch Period RLA/DLA Plot
figure(); hold on;
subplot(2,1,1);
plot(toDateNum(s.X(total_logical)), wrapTo360(rad2deg(s.RLA(total_logical))), 'm-');
dateformat = 2; datetick('x', dateformat); xtickangle(45);
xlim([toDateNum(launch_period_MJDs_launch(1)), toDateNum(launch_period_MJDs_launch(end))]);
ylim([0 180]); 
xlabel('Earth Departure Date'); ylabel('RLA [deg]');
title('RLA throughout Launch Period');
grid on;

subplot(2,1,2); grid on; 
plot(toDateNum(s.X(total_logical)), wrapTo180(rad2deg(s.DLA(total_logical))), 'b-');
dateformat = 2; datetick('x', dateformat); xtickangle(45);
xlim([toDateNum(launch_period_MJDs_launch(1)), toDateNum(launch_period_MJDs_launch(end))]);
patch([700000 800000 800000 700000], [28.57 28.57 90 90], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
patch([700000 800000 800000 700000], [-28.57 -28.57 -90 -90], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
ylim([-45 45]);
xlabel('Earth Departure Date'); ylabel('DLA [deg]');
title('DLA throughout Launch Period');
grid on;

% legend('RLA', 'DLA');
% xlabel('Earth Departure Date'); ylabel('RLA/DLA [deg]');
% title('RLA and DLA throughout Launch Period');
