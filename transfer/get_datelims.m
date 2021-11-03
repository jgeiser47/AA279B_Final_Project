function [x_beg, x_end, y_beg, y_end] = get_datelims()
% Purpose: Get beginning/end MJDs for X and Y axes of search space

% 2033 launch window
x_beg = cal_to_MJD(2032, 8, 1, 0, 0, 0);
x_end = cal_to_MJD(2034, 1, 5, 0, 0, 0);
y_beg = cal_to_MJD(2033, 4, 1, 0, 0, 0);
y_end = cal_to_MJD(2035, 9, 1, 0, 0, 0);
end