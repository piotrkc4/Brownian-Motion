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
%sample size and timestep
dt = 10e-9; %10 nanoseconds
n=100;
%%%%%%%%%%%%%%%%%%%%%%%%
%generating wi and timesteps
wi= randn(n,1);
t = timestep(0,tau,n);%calculating time range comparable with tau
t_over_tau = t./tau; % calculating time constant
%%%%%%%%%%%%%%%%%%%%%%%%
%calculating x (position)
x_inertia = distance_with_inertia(dt, wi, g, m, T, kb, n); %calculating x with inertia
x_inertia = x_inertia/1e-09; %change from meters to nanometers
x_without_inertia = distance_no_inertia(dt, wi, D, n); %calculating x without inertia
x_without_inertia = x_without_inertia/1e-09; %change from meters to nanometers
%%%%%%%%%%%%%%%
%Plots
figure
plot(t_over_tau,x_inertia)
hold on
plot(t_over_tau,x_without_inertia)
ylabel('Xi, (nm)')
xlabel('t/tau')
legend('Brownian movement in inertial regime','Brownian movement in diffusive regime')
%%%%%%
%Calculations of autocorrelation function
Cv_without_inertia = calculate_Cv(dt,x_without_inertia,n);%claclulate Cv without inertia
Cv_with_inertia = calculate_Cv(dt,x_inertia,n); %claclulate Cv with inertia
t3 = timestep(0,6 * tau,length(Cv_with_inertia)); %calculate time range for Cv plot
t_over_tau3 = t3./tau; %calculate time range for Cv plot
%%%%%%%%%%%%%%%%%
%Plots
figure
plot(t_over_tau3, Cv_with_inertia)
hold on
plot(t_over_tau3, Cv_without_inertia)
ylabel('Cv autocorrelation functio')
xlabel('t over tau')
legend('Brownian movement in inertial regime','Brownian movement in diffusive regime')
ylim([0 1]);

