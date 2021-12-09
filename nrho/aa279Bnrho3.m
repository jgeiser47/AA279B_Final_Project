%% AA 279B Final Project
%% Walter Manuel
%% aa279Bnrho3.m
%% 2 November 2021
%% This file plots the northern and southern L2 families

close all
clear
clc
format long
%% Load simulated ICs, define constants

% Saved halo family ICs and period values found in ~nrho2.m
load SMnorthL2.mat

mu_mars = 42828.375816;
mu_sun = 1.32712440018e11;
mu_star = mu_mars/(mu_mars + mu_sun);

%% Plots of Northern L2 Orbits


options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12, 'Events',@ZeroCrossingEvent);

figure;
grid on
hold on
xlabel('x (DU)')
ylabel('y (DU)')
zlabel('z (DU)')
plot3(1-mu_star,0,0,'ro')
plot3(1.004763106761945,0,0,'g*')
% L1 =  0.995251326823620

% L2 =  1.004763106761945

for i = 1:length(t_fam)
    
    [t_out, X_out, ~,~,~] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star),[0 t_fam(i)], X_fam(:,i), options);

    plot3(X_out(:,1),X_out(:,2),X_out(:,3),'b')
    plot3(X_out(:,1),-X_out(:,2),X_out(:,3),'b')

end
legend('Mars','Sun-Mars L2 Point','Halo Orbits')
axis equal
axis square
title('Northern Halo Orbits, Sun-Mars L2 Point')

%% Plots of Southern L2 Orbits

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12, 'Events',@ZeroCrossingEvent);

figure;
grid on
hold on
xlabel('x (DU)')
ylabel('y (DU)')
zlabel('z (DU)')
plot3(1-mu_star,0,0,'ro')
plot3(1.004763106761945,0,0,'g*')
% L1 =  0.995251326823620

% L2 =  1.004763106761945

for i = 1:length(t_fam)
    
    [t_out, X_out, ~,~,~] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star),[0 t_fam(i)], X_fam(:,i), options);

    plot3(X_out(:,1),X_out(:,2),-X_out(:,3),'b')
    plot3(X_out(:,1),-X_out(:,2),-X_out(:,3),'b')

end
legend('Mars','Sun-Mars L2 Point','Halo Orbits')
axis equal
axis square
title('Southern Halo Orbits, Sun-Mars L2 Point')






