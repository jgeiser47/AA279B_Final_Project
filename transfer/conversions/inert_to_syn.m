function r_syn = inert_to_syn(MJDs, r_inert)

% Assumes heliocentric orbit
mu_Sun = 1.3271244004193938e11;

% Reference values
a_Mars = 227953016; 
n_Mars = sqrt(mu_Sun/(a_Mars^3));
AOP_Mars = -1.27848943960252;
M0 = 1.06623161056081;
T0 = cal_to_MJD(2034, 1, 28, 0, 0, 0);

% Rotation matrix to account for slight inclination
RAAN_Mars = 0.863237189337131;
R_0 = [+cos(RAAN_Mars), +sin(RAAN_Mars), 0; ...
       -sin(RAAN_Mars), +cos(RAAN_Mars), 0; ...
                     0,               0, 1];

% Rotation matrix to account for slight inclination
i_Mars = 0.0322348483513284;
R_1 = [1,            0,            0; ...
       0, +cos(i_Mars), +sin(i_Mars); ...
       0, -sin(i_Mars), +cos(i_Mars)];
   
% Setup of variables
N = length(MJDs);
r_syn = zeros(size(r_inert));

for i = 1:N
    MJD = MJDs(i);
    ang = (M0 + AOP_Mars) + n_Mars * ((MJD - T0)*86400);
    R_2 = [+cos(ang), +sin(ang), 0; ...
           -sin(ang), +cos(ang), 0; ...
                   0,         0, 1];
               
    r_syn_i = R_2 * R_1 * R_0 * r_inert(i,1:3)';
    r_syn(i,1:3) = r_syn_i';
end
end