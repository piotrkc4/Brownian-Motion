clc
clear all
close all
dt1= 0.1;
dt2= 0.5;
dt3= 1; %define 3 values of dt
n1 = 200;
n2 = 50;
n3 = 30;
wi1= randn(n1,1);
wi2= randn(n2,1);
wi3= randn(n3,1);   %three sets of random numbers from Gaussian
Wi1 = whitenoise(wi1,dt1);
Wi2 = whitenoise(wi2,dt2);
Wi3 = whitenoise(wi3,dt3); %Three sets of Whitenoise parameters
startValue = 0;
endValue = 30;
t1 = timestep(startValue,endValue,n1); %Evaluating timesteps
t2 = timestep(startValue,endValue,n2);
t3 = timestep(startValue,endValue,n3);
figure
scatter(t1,Wi1,'.')
xlabel('time (s)')
ylabel('Wi')
title('dt=0.1s')
ylim([-8 8]);  %set limits
figure
scatter(t2,Wi2,'.')
xlabel('time (s)')
ylabel('Wi')
title('dt=0.5s')
ylim([-8 8]); %set limits
figure
scatter(t3,Wi3,'.')
xlabel('time (s)')
ylabel('Wi')
title('dt=1s')
ylim([-8 8]); %set limits
x1= distance(n1, wi1, dt1); %calculate x with time
x2= distance(n2, wi2, dt2);
x3= distance(n3, wi3, dt3);
figure
plot(t1,x1);
xlabel('time (s)')
ylabel('xi')
title('dt=0.1s')
figure
plot(t2,x2);
xlabel('time (s)')
ylabel('xi')
title('dt=0.5s')
figure
plot(t3,x3);
xlabel('time (s)')
ylabel('xi')
title('dt=1s')
Mean = [mean(x1) mean(x2) mean(x3)]; %calculate means
Var1 = averaging(x1,Mean(:,1));
Var2 = averaging(x2,Mean(:,2));
Var3 = averaging(x3,Mean(:,3));









