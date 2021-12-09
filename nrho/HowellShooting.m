function [controls] = HowellShooting(phi,X,mu,fix0)
%Single Shooting method for perpendicular crossings

x = X(1);
y = X(2);
z = X(3);
ydot = X(5);

d1 = sqrt((x + mu)^2 + y^2 + z^2);
d2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

dU_dx = x + ((mu + x)*(mu - 1))/(d1^3) - (mu*(mu + x - 1))/(d2^3);
%dU_dy = y - (mu*y)/(d2^3) + (y*(mu - 1))/(d1^3);
dU_dz = (z*(mu - 1))/(d1^3) - (mu*z)/(d2^3);

xddot = 2*X(5) + dU_dx;
%Xp(5) = -2*X(4) + dU_dy; % y-doubledot
zddot = dU_dz;

errors = [0; 0] - [X(4); X(6)];

switch fix0
    case 1 % keep x0 fixed
        sub_phi = [phi(4,3) phi(4,5);
                   phi(6,3) phi(6,5)];
               
        temp = (1/ydot)*[xddot; zddot]*[phi(2,3) phi(2,5)];
        M = sub_phi - temp;
    case 3 % keep z0 fixed
        sub_phi = [phi(4,1) phi(4,5);
                   phi(6,1) phi(6,5)];
               
        temp = (1/ydot)*[xddot; zddot]*[phi(2,1) phi(2,5)];
        M = sub_phi - temp;
end

controls = M\errors;

end








