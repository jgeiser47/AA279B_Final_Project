function [x_beg, x_end, y_beg, y_end] = get_datelims()
% Purpose: Get beginning/end MJDs for X and Y axes of search space

% 2020-2021 launch window (Mars 2020 dates: July 30, 2020 - Feb 18, 2021)
x_beg = cal_to_MJD(2020, 1, 1, 0, 0, 0);
x_end = cal_to_MJD(2021, 1, 5, 0, 0, 0);
y_beg = cal_to_MJD(2020, 8, 1, 0, 0, 0);
y_end = cal_to_MJD(2022, 9, 1, 0, 0, 0);
end