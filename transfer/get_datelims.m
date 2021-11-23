function [x_beg, x_end, y_beg, y_end] = get_datelims()
% Purpose: Get beginning/end MJDs for X and Y axes of search space

% 2033 launch window
x_beg = cal_to_MJD(2032, 12, 1, 0, 0, 0);
x_end = cal_to_MJD(2033, 8, 14, 0, 0, 0);
y_beg = cal_to_MJD(2033, 7, 1, 0, 0, 0);
y_end = cal_to_MJD(2034, 7, 20, 0, 0, 0);
end