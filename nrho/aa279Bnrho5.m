%% AA 279B Final Project
%% Walter Manuel
%% aa279Bnrho5.m
%% 29 November 2021
%% This file is mainly for plotting

close all
clear
clc
format long
%% Load simulated ICs, define constants

s1 = load('SMnorthL2b.mat'); %100 northern halo orbits
t_fam = s1.t_fam;
X_fam = s1.X_fam;
s2 = load('SMnorthL2c.mat'); % 10 more nrho orbits
t_fam2 = s2.t_fam2;
X_fam2 = s2.X_fam2;

load('HaloEigenvalues.mat');

mu_mars = 42828.375816;
mu_sun = 1.32712440018e11;
mu_sun = 1.3271244004193938e11;
mu_star = mu_mars/(mu_mars + mu_sun);

Lstar = 227953016;
Tstar = sqrt(Lstar^3/(mu_mars + mu_sun));

%% Stability Check

Evalsort = sort(Eval,2,'ComparisonMethod','real');
v1 = [];
v2 = [];

% had to check these manually since Matlab keeps changing the eigenvalue
% order
for i = 1:length(Eval)
    if i < 11
        v1(i) = (Evalsort(i,1) + Evalsort(i,6))/2; % real (mostly)
        v2(i) = (Evalsort(i,2) + Evalsort(i,3))/2; % complex
    elseif i < 45
        v1(i) = (Evalsort(i,3) + Evalsort(i,6))/2; % positive
        v2(i) = (Evalsort(i,1) + Evalsort(i,2))/2; % negative
    elseif i < 100
        v1(i) = (Evalsort(i,3) + Evalsort(i,4))/2;
        v2(i) = (Evalsort(i,1) + Evalsort(i,2))/2;
    elseif i < 111
        v1(i) = (Evalsort(i,4) + Evalsort(i,5))/2;
        v2(i) = (Evalsort(i,1) + Evalsort(i,2))/2;
    end
end

%% Stability plots
figure;
plot(([X_fam(1,:),X_fam2(1,:)]-1)*Lstar,abs(v1))
hold on
plot(([X_fam(1,:),X_fam2(1,:)]-1)*Lstar,abs(v2))
plot([7.04e5 7.04e5],[0,600],'--')
grid on
legend('Stability Index 1','Stability Index 2','NRHO Transition Point')
xlabel('Distance past Mars along x-axis in synodic frame (km)')
ylabel('Stability Index (absolute value)')
title('Stability Indices for Northern L2 Sun-Mars Halo Orbits')

figure;
plot(abs(v1))
hold on
plot(abs(v2))
grid on
% 44 = nrho start
%% Plots of Northern L2 NRHOs

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
    
    if i < 44
    plot3(X_out(:,1),X_out(:,2),X_out(:,3),'b')
    plot3(X_out(:,1),-X_out(:,2),X_out(:,3),'b')
    else 
        plot3(X_out(:,1),X_out(:,2),X_out(:,3),'m')
        plot3(X_out(:,1),-X_out(:,2),X_out(:,3),'m')
    end

end
legend('Mars','Sun-Mars L2 Point','Normal Halo Orbits')
axis equal
axis square
title('Northern Halo Orbits with NRHOs Highlighted')

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

for i = 44:length(t_fam)
    
    [t_out, X_out, ~,~,~] = ode113(@(t,x)CR3BP_EOM_3(t,x,mu_star),[0 t_fam(i)], X_fam(:,i), options);

        plot3(X_out(:,1),X_out(:,2),X_out(:,3),'m')
        plot3(X_out(:,1),-X_out(:,2),X_out(:,3),'m')


end
legend('Mars','Sun-Mars L2 Point','NRHOs')
axis equal
axis square
title('Northern L2 NRHO Family')

%% Nominal NRHO plot

s3 = load('nominal_north_NRHO.mat');
X_nrho = s3.X_nrho;
t_nrho = s3.t_nrho;

figure;
grid on
hold on
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
plot3((1-mu_star)*Lstar,0,0,'ro')
plot3((1.004763106761945)*Lstar,0,0,'g*')

plot3(X_nrho(:,1)*Lstar,X_nrho(:,2)*Lstar,X_nrho(:,3)*Lstar,'m')
legend('Mars','Sun-Mars L2 Point','Nominal NRHO')
axis equal
axis square
title('Nominal NRHO in Synodic Frame')


%% Nominal NRHO properties

period = (t_nrho(end)*Tstar)
period = period/60/60/24

apoapsis = X_nrho(1,3)*Lstar
periapsis = X_nrho(207,3)*Lstar


v1(82)
v2(82)

%% Resonance
Tmars = 2*pi;
resonance = Tmars./[t_fam t_fam2]/2

figure;
plot(([X_fam(1,:),X_fam2(1,:)]-1)*Lstar,resonance)
grid on
title('Resonance Ratios of Sun-Mars L2 Halo Orbits')
xlabel('Distance past Mars along x-axis in synodic frame (km)')
ylabel('Revolutions of Halo Orbit per Martian Synodic Period')



