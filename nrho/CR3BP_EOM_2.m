function [Xp] = CR3BP_EOM_2(t,X,mu1,mu2,d)
% X = x,y,z,x-dot,y-dot,z-dot
% no state transition matrix
% Inertial Frame
% Dimensional (km,s)

x = X(1);
y = X(2);
z = X(3);

r3 = [x; y; z];

mu = mu2/(mu1 + mu2);
ws = sqrt((mu1 + mu2)/(d^3));

R1 = mu*d;
R2 = (1-mu)*d;

r1 = [-R1*cos(ws*t); -R1*sin(ws*t); 0];
r2 = [R2*cos(ws*t); R2*sin(ws*t); 0];

r13 = r3 - r1;
r23 = r3 - r2;

r13x = r13(1);
r13y = r13(2);
r13z = r13(3);

r23x = r23(1);
r23y = r23(2);
r23z = r23(3);

Xp = ones(6,1); % [x1 x2 x3 x1-dot x2-dot x3-dot];
Xp(1) = X(4); % x-dot;
Xp(2) = X(5); % y-dot;
Xp(3) = X(6); % z-dot;

Xp(4) = -mu1*r13x/(norm(r13)^3) - mu2*r23x/(norm(r23)^3); % x-doubledot
Xp(5) = -mu1*r13y/(norm(r13)^3) - mu2*r23y/(norm(r23)^3); % y-doubledot
Xp(6) = -mu1*r13z/(norm(r13)^3) - mu2*r23z/(norm(r23)^3); % z-doubledot

end



