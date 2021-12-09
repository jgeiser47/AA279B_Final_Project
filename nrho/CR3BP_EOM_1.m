function [Xp] = CR3BP_EOM_1(t,X,mu)
% X = x,y,z,x-dot,y-dot,z-dot
% no state transition matrix
% Synodic Frame
% Non-Dimensional

x = X(1);
y = X(2);
z = X(3);

d1 = sqrt((x + mu)^2 + y^2 + z^2);
d2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
% U = 0.5*(x^2 + y^2) + (1 - mu)/d1 + mu/d2;

dU_dx = x + ((mu + x)*(mu - 1))/(d1^3) - (mu*(mu + x - 1))/(d2^3);
dU_dy = y - (mu*y)/(d2^3) + (y*(mu - 1))/(d1^3);
dU_dz = (z*(mu - 1))/(d1^3) - (mu*z)/(d2^3);

Xp = ones(6,1); % [x1 x2 x3 x1-dot x2-dot x3-dot];
Xp(1) = X(4); % x-dot;
Xp(2) = X(5); % y-dot;
Xp(3) = X(6); % z-dot;

Xp(4) = 2*X(5) + dU_dx; % x-doubledot
Xp(5) = -2*X(4) + dU_dy; % y-doubledot
Xp(6) = dU_dz; % z-doubledot

end



