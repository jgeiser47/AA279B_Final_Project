function fx = f(x1, x2)
% Purpose: Wrapper for optimization problem. Returns f(x) given inputs x1
%          and x2, where x1 is the Earth departure MJD and x2 is the Mars 
%          arrival MJD

% Calculate total dV
data = solve_lambert(x1, x2);
fx = data.C3_1;
end