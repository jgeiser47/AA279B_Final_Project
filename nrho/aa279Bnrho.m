%% AA 279B Final Project
%% Walter Manuel
%% aa279Bnrho.m
%% 2 November 2021
%% This file finds a foundational L2 Halo Orbit

close all
clear
clc
format long
%% Original Halo Orbit

mu_mars = 42828.375816;
%mu_sun = 1.32712440018e11;
mu_sun = 1.3271244004193938e11;
mu_star = mu_mars/(mu_mars + mu_sun);

% this IC taken from cited paper
% Sun-Mars northern L2
X0 = [1.004291664675 + 0.001; 0; -.00043806867; 0; -.002464489308; 0];

tspan = linspace(0,5,1000);

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12, 'Events',@ZeroCrossingEvent);
[t_out, X_out, tE,yE,iE] = ode113(@(t,x)CR3BP_EOM_1(t,x,mu_star), tspan, X0, options);

figure
plot3(X_out(:,1),X_out(:,2),X_out(:,3))
grid on
xlabel('x')
ylabel('y')
zlabel('z')

%% Use shooting method to correct original Halo
phi0 = eye(6);
phi0 = reshape(phi0,36,1);
X_0 = [X0; phi0];
del_x = 0;
del_z = 0;
del_ydot = 0;
count = 1;
err = 1;

figure;
while count <= 100 && err > 1e-8
    count = count + 1;
    
    X_0 = X_0 + [del_x; 0; del_z; 0; del_ydot; 0; zeros(36,1)];

    [t_out, X_out, tE,yE,iE] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star), tspan, X_0, options);
    
    phi = X_out(:,7:42);
    phi = reshape(phi,[],6,6);
    phi = permute(phi,[2 3 1]);

    phi_end = phi(:,:,end);
    X_end = X_out(end,1:6);
     
    controls = HowellShooting(phi_end,X_end,mu_star,1);
    
    del_z = controls(1);
    del_ydot = controls(2);
    targets = [0 - X_end(4); 0 - X_end(6)];
    err = norm(targets);
    

    plot3(X_out(:,1),X_out(:,2),X_out(:,3))
    grid on
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z')
end

%% Plot corrected foundational Halo orbit

[t_out, X_out] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star), tspan, X_0);


figure;
plot3(X_out(:,1),X_out(:,2),X_out(:,3))
hold on
grid on
axis equal
%plot3(X_out(:,1),-X_out(:,2),X_out(:,3))

% L1 =  0.995251326823620

% L2 =  1.004763106761945


















