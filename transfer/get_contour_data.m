function s = get_contour_data()
% Purpose: Get data for contour lines by running Lambert's algorithm across
%          search space. This function helps prevent running this over and
%          over again for each contour plot. 

% Output struct
s = struct(); 

% Get axis limits
[x_beg, x_end, y_beg, y_end] = get_datelims();

% How fine the mesh grid is
grid_dt = 1; 
N = 100;

% Get arrays of MJDs
x = [x_beg : grid_dt : x_end]; % x = linspace(x_beg, x_end, N); 
y = [y_beg : grid_dt : y_end]; % y = linspace(y_beg, y_end, N); 

% Get mesh grid and setup Z (aka f(x)) variable
[X,Y] = meshgrid(x,y);
Z = zeros(size(X));
RLA = zeros(size(X));
DLA = zeros(size(X));
v1_out = zeros([size(X), 3]);

% Populate f(x) values given inputs x1 and x2
for i = 1:size(X,1)
    for j = 1:size(X,2)
%         data_ij = solve_lambert(X(i,j), Y(i,j));
%         Z(i,j) = data_ij.C3_1;
%         data(i,j) = data_ij;

        data_ij = solve_lambert(X(i,j), Y(i,j));
        Z(i,j) = data_ij.C3_1;
        RLA(i,j) = data_ij.RLA;
        DLA(i,j) = data_ij.DLA;
        v1_out(i,j,:) = data_ij.v1_out;
    end
end

% Constraint lines (0 < TOF < 2 years)
c1_x = linspace(x_beg, x_end+100, N);
c1_y = c1_x;
c2_y = c1_x + 365*2;

% Add values to struct
s.N = N;
s.x_beg = x_beg; s.x_end = x_end; 
s.y_beg = y_beg; s.y_end = y_end;
s.X = X; s.Y = Y; s.Z = Z;
s.c1_x = c1_x; s.c1_y = c1_y; s.c2_y = c2_y;
% s.levels = [5,5.5,6,6.5,7,7.5,8,8.5,9,10,11,12,13,14,15];
s.levels = [7,8,9,10,12,14,16,18,20,25,30];
s.RLA = RLA; s.DLA = DLA; s.v1_out = v1_out;
end