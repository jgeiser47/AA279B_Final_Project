function [Xp] = CR3BP_EOM_3(t,X,mu)
% X = x,y,z,x-dot,y-dot,z-dot,phi
% includes state transition matrix
% Synodic Frame
% Non-Dimensional

x = X(1);
y = X(2);
z = X(3);

phi = X(7:42); % Transition Matrix
phi = reshape(phi,6,6);
phi_dot = ones(6);

d1 = sqrt((x + mu)^2 + y^2 + z^2);
d2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
% U = 0.5*(x^2 + y^2) + (1 - mu)/d1 + mu/d2;

dU_dx = x + ((mu + x)*(mu - 1))/(d1^3) - (mu*(mu + x - 1))/(d2^3);
dU_dy = y - (mu*y)/(d2^3) + (y*(mu - 1))/(d1^3);
dU_dz = (z*(mu - 1))/(d1^3) - (mu*z)/(d2^3);


d2U_dx2 = (mu - 1)/(d1^3) - mu/(d2^3) + (3*mu*(mu + x - 1)^2)/(d2^5) - (3*(mu + x)^2*(mu - 1))/(d1^5) + 1;
d2U_dxdy = (3*mu*y*(mu + x - 1))/(d2^5) - (3*y*(mu + x)*(mu - 1))/(d1^5);
d2U_dxdz = (3*mu*z*(mu + x - 1))/(d2^5) - (3*z*(mu + x)*(mu - 1))/(d1^5);

d2U_dy2 = (mu - 1)/(d1^3) - mu/(d2^3) - (3*y^2*(mu - 1))/(d1^5) + (3*mu*y^2)/(d2^5) + 1;
d2U_dydz = (3*mu*y*z)/(d2^5) - (3*y*z*(mu - 1))/(d1^5);

d2U_dz2 = (mu - 1)/(d1^3) - mu/(d2^3) - (3*z^2*(mu - 1))/(d1^5) + (3*mu*z^2)/(d2^5);

mat_Z = zeros(3);
mat_I = eye(3);
mat_W = 2 * [0 1 0; -1 0 0; 0 0 0];
mat_U = [d2U_dx2    d2U_dxdy   d2U_dxdz; 
         d2U_dxdy   d2U_dy2    d2U_dydz; 
         d2U_dxdz   d2U_dydz   d2U_dz2];

F = [mat_Z   mat_I; 
     mat_U   mat_W];

Xp = ones(6,1); % [x1 x2 x3 x1-dot x2-dot x3-dot];
Xp(1) = X(4); % x-dot;
Xp(2) = X(5); % y-dot;
Xp(3) = X(6); % z-dot;

Xp(4) = 2*X(5) + dU_dx; % x-doubledot
Xp(5) = -2*X(4) + dU_dy; % y-doubledot
Xp(6) = dU_dz; % z-doubledot

phi_dot = F*phi;
phi_dot = reshape(phi_dot,36,1);

Xp = [Xp; phi_dot];

end



