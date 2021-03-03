clc
clear all
close all
%Data entry
R = 1e-06; %in meters
m = 1.1e-14; %in kilograms
V = 0.001; % in Ns/m2
g = 6 * pi * V * R; % Ns/m2 * m = Ns/m
T = 300; %in K
tau = 6e-7; % in seconds
kb = 1.38e-23; %in m2*kg/s2*K
D = (kb*T)/g; %diffusion coeff m2*kg/s2 *  1/kg = m2/s2
%%%%%%%%%%%%%%%%%%%%%%%%%
%Gaussian generation
dt = 10e-9; %10 nanoseconds
n=1000;
%%%%%%%%%%%%%%%%%%%%%%%%
%generating wi and timesteps
wi= randn(n,1);
startValue = 0;
endValue = 100 * tau;
nElements = n;
t = timestep(startValue,endValue,nElements);%calculating time range comparable with tau
t_over_tau = t./tau; % calculating time constant
t2 = timestep(0,100 * tau,316);
t_over_tau2 = t2./tau;
%%%%%%%%%%%%%%%%%%%%%%%%%
%calculating x
x_inertia = distance_with_inertia(dt, wi, g, m, T, kb, n); %calculating x with inertia
x_inertia = x_inertia/1e-09; %change from meters to nanometers
x_without_inertia = distance_no_inertia(dt, wi, D, n); %calculating x without inertia
x_without_inertia = x_without_inertia/1e-09; %change from meters to nanometers
%%%%%%%%%%%%%%%%
%Plots
figure
plot(t_over_tau,x_inertia)
ylabel('Xi, (nm)')
xlabel('t/tau')
ylim([-5 5]); %ylimits
xlim([0 100]); %x limits
hold on
plot(t_over_tau,x_without_inertia)
ylabel('Xi, (nm)')
xlabel('t/tau')
legend('Brownian movement in inertial regime','Brownian movement in diffusive regime')
%%%%%%%%%%%%%%%%%%
%Calculate MSD
MSD_without_inertia = mean_square_displacement(x_without_inertia,n); %calculate mean square displacement for x_without_inertia
MSD_inertia = mean_square_displacement(x_inertia,n);
figure
loglog(t_over_tau2,MSD_inertia)
hold on
loglog(t_over_tau2,MSD_without_inertia)
xlim([0 100])