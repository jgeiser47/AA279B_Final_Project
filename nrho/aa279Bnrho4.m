%% AA 279B Final Project
%% Walter Manuel
%% aa279Bnrho3.m
%% 2 November 2021
%% This file does the bulk of the work to find NRHOs from the Halo family

close all
clear
clc
format long
%% Load simulated ICs, define constants

load SMnorthL2.mat

mu_mars = 42828.375816;
%mu_sun = 1.32712440018e11;
mu_sun = 1.3271244004193938e11;
mu_star = mu_mars/(mu_mars + mu_sun);

Lstar = 227953016;
Tstar = sqrt(Lstar^3/(mu_mars + mu_sun));

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
impactCheck = [];

for i = 1:length(t_fam)
    
    [t_out, X_out, ~,~,~] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star),[0 t_fam(i)], X_fam(:,i), options);

    plot3(X_out(:,1),X_out(:,2),X_out(:,3),'b')
    plot3(X_out(:,1),-X_out(:,2),X_out(:,3),'b')
    
    impactCheck(i,1) = X_out(end,3)*Lstar;

end
legend('Mars','Sun-Mars L2 Point','Halo Orbits')
axis equal
axis square
title('Northern Halo Orbits, Sun-Mars L2 Point')

%% Plots of Southern L2 Orbits

options1 = odeset('RelTol', 1e-9, 'AbsTol', 1e-12, 'Events',@ZeroCrossingEvent);
options2 = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
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

for i = length(t_fam)
    
    [t_out, X_out] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star),[0 t_fam(i)], X_fam(:,i), options2);

    plot3(X_out(:,1),X_out(:,2),-X_out(:,3),'b')
    %plot3(X_out(:,1),-X_out(:,2),-X_out(:,3),'b')

end
legend('Mars','Sun-Mars L2 Point','Halo Orbits')
axis equal
axis square
title('Southern Halo Orbits, Sun-Mars L2 Point')

%% First attempt to calculate stability indices and determine NRHOs
Eval = [];
for i = 1:length(t_fam)
    
    [t_out, X_out] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star),[0 t_fam(i)]*2, X_fam(:,i), options2);

     phi_out = X_out(:,7:42);
     phi_out = reshape(phi_out,[],6,6);
     phi_out = permute(phi_out,[2 3 1]);
     monodromy = phi_out(:,:,end);
     Eval(i,:) = eig(monodromy);
     v(i,:) = CR3BP_Stability(monodromy);

end

figure;
plot(v(:,1))
hold on
plot(v(:,2))
grid on

%% Generated 10 more Halo Orbits that weren't that useful
X_0 = X_fam(:,end);
del_x = 0;
del_z = 0;
del_ydot = 0;

err = 1;

tspan = linspace(0,5,5000);

nHalos = 10;
dx = -0.0002/nHalos;
figure;

for i = 1:nHalos
    X_0 = X_0 + [dx; 0; 0; 0; 0; 0; zeros(36,1)];
    count = 1;
    err = 1;
    while count <= 100 && err > 1e-8
        count = count + 1;

        X_0 = X_0 + [del_x; 0; del_z; 0; del_ydot; 0; zeros(36,1)];

        [t_out, X_out, tE,yE,iE] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star), tspan, X_0, options1);

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
    t_fam2(i) = t_out(end);
    X_fam2(:,i) = X_0;

    impactCheck2(i,1) = X_out(end,3)*Lstar;
    
    plot3(X_out(:,1),X_out(:,2),X_out(:,3))
    grid on
    hold on
    axis equal
    axis square
    xlabel('x')
    ylabel('y')
    zlabel('z')

end

%% Attempt to calculate stability indices and determine NRHOs from new 10 orbits
for i = 1:length(t_fam2)
    
    [t_out, X_out] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star),[0 t_fam2(i)]*2, X_fam2(:,i), options2);

     phi_out = X_out(:,7:42);
     phi_out = reshape(phi_out,[],6,6);
     phi_out = permute(phi_out,[2 3 1]);
     monodromy = phi_out(:,:,end);
     [EvecA,EvalA(i,:)] = eig(monodromy, 'vector');
     v2(100+i,:) = CR3BP_Stability(monodromy);

end
%% Tried to plot stability values to verify them
nrhos = 1:100;
apoapsis = [X_fam(3,nrhos) X_fam2(3,:)]

figure;
plot(v(:,1))
hold on
plot(v(:,2))
grid on

figure;
plot(impactCheck,[v(nrhos,1)])
hold on
plot(impactCheck,[v(nrhos,2)])
grid on

figure;
plot([v(nrhos,1); v2(:,1)])
hold on
plot([v(nrhos,2); v2(:,2)])
grid on

figure;
plot(apoapsis,[v(nrhos,1); v2(:,1)])
hold on
plot(apoapsis,[v(nrhos,2); v2(:,2)])
grid on

%% Plots of NRHO orbits

figure;
grid on
hold on
xlabel('x (DU)')
ylabel('y (DU)')
zlabel('z (DU)')
plot3(1-mu_star,0,0,'ro')
plot3(1.004763106761945,0,0,'g*')
for i = 40:length(t_fam)
    
    [t_out, X_out] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star),[0 t_fam(i)], X_fam(:,i), options2);
    plot3(X_out(:,1),X_out(:,2),-X_out(:,3),'b')
    %plot3(X_out(:,1),-X_out(:,2),-X_out(:,3),'b')
end

for i = 1:length(t_fam2)
    
    [t_out, X_out] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star),[0 t_fam2(i)], X_fam2(:,i), options2);
    plot3(X_out(:,1),X_out(:,2),-X_out(:,3),'b')
    %plot3(X_out(:,1),-X_out(:,2),-X_out(:,3),'b')
end

legend('Mars','Sun-Mars L2 Point','Halo Orbits')
axis equal
axis square
title('Southern Halo Orbits, Sun-Mars L2 Point')

%% Nominal NRHO and plot

X_nom = X_fam(:,82);
t_nom = t_fam(82);
[t_out, X_out] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star),[0 t_nom]*2, X_nom, options2);
    
figure;
plot3(X_out(:,1),X_out(:,2),X_out(:,3),'b')
grid on
hold on
%plot3(X_out(:,1),-X_out(:,2),X_out(:,3),'b')
xlabel('x (DU)')
ylabel('y (DU)')
zlabel('z (DU)')
plot3(1-mu_star,0,0,'ro')
plot3(1.004763106761945,0,0,'g*')
axis equal

% this output went to nominal_north_NRHO.mat
X_nrho = X_out;
t_nrho = t_out;

%% Dimensionalized plots and ICs
JDstart = 2464008.50000;
[r_mars,v_mars] = helioPosVel(4,JDstart);
t0 = atan2(r_mars(2),r_mars(1));
t0 = atan(.9);
for i = 1:length(X_nrho)
    [X_nrhoD(i,:),t_nrhoD(i)] = CR3BPtoHCI(X_nrho(i,1:6),t_nrho(i),t0);
    %norm(X_nrhoD(i,4:6))
end

