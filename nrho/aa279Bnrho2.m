%% AA 279B Final Project
%% Walter Manuel
%% aa279Bnrho2.m
%% 2 November 2021
%% This file generates the family of Halo Orbits

close all
clear
clc
format long
%% Foundational Halo Orbit (found in ~nrho.m)

mu_mars = 42828.375816;
%mu_sun = 1.32712440018e11;
mu_sun = 1.3271244004193938e11;
mu_star = mu_mars/(mu_mars + mu_sun);

%X0 = [1.00529166467500;0;0.00120385999464693;0;-0.00465015945829012;0];
X0 = [1.00529166467500;0;0.001203859989615;0;-0.004650159454831;0];

tspan = linspace(0,5,5000);

% event option will make it stop at the half period point
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12, 'Events',@ZeroCrossingEvent);
[t_out, X_out, tE,yE,iE] = ode113(@(t,x)CR3BP_EOM_1(t,x,mu_star), tspan, X0, options);

% plot foundational halo orbit
figure
plot3(X_out(:,1),X_out(:,2),X_out(:,3))
grid on
xlabel('x')
ylabel('y')
zlabel('z')

%% Use continuation scheme and shooting method to generate Halo family
phi0 = eye(6);
phi0 = reshape(phi0,36,1);
X_0 = [X0; phi0];
del_x = 0;
del_z = 0;
del_ydot = 0;

err = 1;


nHalos = 100;
dx = -0.005/nHalos; % increment towards Mars away from L2 point
figure;

for i = 1:nHalos % increment
    X_0 = X_0 + [dx; 0; 0; 0; 0; 0; zeros(36,1)];
    count = 1;
    err = 1;
    while count <= 100 && err > 1e-8 % shooting method
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



    end
    t_fam(i) = t_out(end);
    X_fam(:,i) = X_0;
    
    % plot each newly found family member
    plot3(X_out(:,1),X_out(:,2),X_out(:,3))
    grid on
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z')
end

%% plot last orbit calculated

%[t_out, X_out] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star), tspan, X_0);


figure;
plot3(X_out(:,1),X_out(:,2),X_out(:,3))
hold on
grid on
axis equal
%plot3(X_out(:,1),-X_out(:,2),X_out(:,3))

% L1 =  0.995251326823620

% L2 =  1.004763106761945


















