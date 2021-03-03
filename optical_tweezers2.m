clc
clear all
close all
%Data entry
R = 1e-06; %in meters
V = 0.001; % in Ns/m2
g = 6 * pi * V * R; % Ns/m2 * m = Ns/m
T = 300; %in K
kb = 1.38e-23; %in m2*kg/s2*K
D = (kb*T)/g; %diffusion coeff m2*kg/s2 *  1/kg = m2/s2
kx = 1e-6; %in N/m
ky = 1e-6; %in N/m
kz = 0.2e-6; %in N/m
k = [kx ky kz];
%%%%%%%%%%%%%%%%%%%%%%%%%
%sample size and step size
phi = g/kx;
dt = 0.001; 
n=1000;
%%%%%%%%%%%%%%%%%%%%%%%%
%generating wi and timesteps
wi= randn(3,n);
t = timestep(0,1,n);%calculating time range comparable with tau
%%%%%%%%%%%%%%%%%%%%%%%%%
%calculating r (position vector)
r_xyz = tweezer_position(dt,k, wi, D, g, n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dependence on Kxy
k_xy = [0.1 0.2 0.5 1 2 3 4 5 6 7 8 9 10]' * 1.0e-6; %N/m
wi_new = [wi(1,:); wi(2,:)];
for j = 1:length(k_xy)
[r_x(j,:),r_y(j,:)] = tweezer_position_XYplane(dt,k_xy(j,:), wi_new, D, g, n);
 r_x(j,:) = r_x(j,:)/1e-9;
 r_y(j,:) = r_y(j,:)/1e-9;
sigma_x(:,j) = var(r_x(j,:));
sigma_y(:,j) = var(r_y(j,:));
Sigma(:,j) = (sigma_x(:,j) + sigma_y(:,j))/2;
end
figure
plot(k_xy,Sigma)
ylabel('sigma_x_y (nm^2)')
xlabel('k_x_y')
%%%%%%%%%%%%%%%%%%%%%
%Prob distributions
%Probability Distribution 1
%x-y plane
%Figure 1 
figure
r_xy = [r_x(1,:);r_y(1,:)];
plot3(r_xy(1,:), r_xy(2,:), zeros(size(r_xy(1:2,:))), 'k.', 'MarkerSize', 1);
xlabel('xt (nm)')
ylabel('yt (nm)')
zlabel('p')
hold on
M1 = mean(r_xy(1,:));
M2 = mean(r_xy(2,:));
S1 = std(r_xy(1,:));
S2 = std(r_xy(2,:));
S = [S1 S2];
M = [M1 M2];
axisMin = M - 4 * S;
axisMax = M + 4 * S;
%elipse
data_zeroMean = bsxfun(@minus, r_xy(1:2,:)', M);
[V,Dm] = eig(data_zeroMean' * data_zeroMean / (size(data_zeroMean, 1)));
[Dm, order] = sort(diag(Dm), 'descend');
Dm = diag(Dm);
V = V(:, order);
V = V * sqrt(Dm);
t = linspace(0, 2 * pi);
e = bsxfun(@plus, 2*V * [cos(t); sin(t)], M');
plot3(...
    e(1,:), e(2,:), ...
    zeros(1, length(e)), 'g-', 'LineWidth', 2);
%Code for histagram
maxP = 0;
for side = 1:2
    % Calculate the histogram.
    p = [0 hist(r_xy(side,:), 20) 0];
    p = p / sum(p);
    maxP = max([maxP p]);
    dx = (axisMax(side) - axisMin(side)) / numel(p) / 2.3;
    p2 = [zeros(1,numel(p)); p; p; zeros(1,numel(p))]; p2 = p2(:);
    x = linspace(axisMin(side), axisMax(side), numel(p));
    x2 = [x-dx; x-dx; x+dx; x+dx]; x2 = max(min(x2(:), axisMax(side)), axisMin(side));

    % Calculate the curve.
    nPtsCurve = numel(p) * 10;
    xx = linspace(axisMin(side), axisMax(side), nPtsCurve);

    % Plot the curve and the histogram.
    if side == 1
        plot3(xx, ones(1, nPtsCurve) * axisMax(3 - side), spline(x,p,xx), 'r-', 'LineWidth', 2);
        plot3(x2, ones(numel(p2), 1) * axisMax(3 - side), p2, 'k-', 'LineWidth', 1);
    else
        plot3(ones(1, nPtsCurve) * axisMax(3 - side), xx, spline(x,p,xx), 'b-', 'LineWidth', 2);
        plot3(ones(numel(p2), 1) * axisMax(3 - side), x2, p2, 'k-', 'LineWidth', 1);
    end

end
grid on
hold off
% for k =1fN
r_xy = [r_x(4,:);r_y(4,:)];
figure
plot3(r_xy(1,:), r_xy(2,:), zeros(size(r_xy(1:2,:))), 'k.', 'MarkerSize', 1);
xlabel('xt (nm)')
ylabel('yt (nm)')
zlabel('p')
hold on
M1 = mean(r_xy(1,:));
M2 = mean(r_xy(2,:));
S1 = std(r_x(1,:));
S2 = std(r_xy(2,:));
S = [S1 S2];
M = [M1 M2];
axisMin = M - 4 * S;
axisMax = M + 4 * S;
%elipse
data_zeroMean = bsxfun(@minus, r_xy(1:2,:)', M);
[V,Dm] = eig(data_zeroMean' * data_zeroMean / (size(data_zeroMean, 1)));
[Dm, order] = sort(diag(Dm), 'descend');
Dm = diag(Dm);
V = V(:, order);
V = V * sqrt(Dm);
t = linspace(0, 2 * pi);
e = bsxfun(@plus, 2*V * [cos(t); sin(t)], M');
plot3(...
    e(1,:), e(2,:), ...
    zeros(1, length(e)), 'g-', 'LineWidth', 2);
%Code for histagram
maxP = 0;
for side = 1:2
    % Calculate the histogram.
    p = [0 hist(r_xy(side,:), 20) 0];
    p = p / sum(p);
    maxP = max([maxP p]);
    dx = (axisMax(side) - axisMin(side)) / numel(p) / 2.3;
    p2 = [zeros(1,numel(p)); p; p; zeros(1,numel(p))]; p2 = p2(:);
    x = linspace(axisMin(side), axisMax(side), numel(p));
    x2 = [x-dx; x-dx; x+dx; x+dx]; x2 = max(min(x2(:), axisMax(side)), axisMin(side));

    % Calculate the curve.
    nPtsCurve = numel(p) * 10;
    xx = linspace(axisMin(side), axisMax(side), nPtsCurve);

    % Plot the curve and the histogram.
    if side == 1
        plot3(xx, ones(1, nPtsCurve) * axisMax(3 - side), spline(x,p,xx), 'r-', 'LineWidth', 2);
        plot3(x2, ones(numel(p2), 1) * axisMax(3 - side), p2, 'k-', 'LineWidth', 1);
    else
        plot3(ones(1, nPtsCurve) * axisMax(3 - side), xx, spline(x,p,xx), 'b-', 'LineWidth', 2);
        plot3(ones(numel(p2), 1) * axisMax(3 - side), x2, p2, 'k-', 'LineWidth', 1);
    end

end
grid on
hold off
%for k = 5fN
r_xy = [r_x(8,:);r_y(8,:)];
figure
plot3(r_xy(1,:), r_xy(2,:), zeros(size(r_xy(1:2,:))), 'k.', 'MarkerSize', 1);
xlabel('xt (nm)')
ylabel('yt (nm)')
zlabel('p')
hold on
M1 = mean(r_xy(1,:));
M2 = mean(r_xy(2,:));
S1 = std(r_xy(1,:));
S2 = std(r_xy(1,:));
S = [S1 S2];
M = [M1 M2];
axisMin = M - 4 * S;
axisMax = M + 4 * S;
%elipse
data_zeroMean = bsxfun(@minus, r_xy(1:2,:)', M);
[V,Dm] = eig(data_zeroMean' * data_zeroMean / (size(data_zeroMean, 1)));
[Dm, order] = sort(diag(Dm), 'descend');
Dm = diag(Dm);
V = V(:, order);
V = V * sqrt(Dm);
t = linspace(0, 2 * pi);
e = bsxfun(@plus, 2*V * [cos(t); sin(t)], M');
plot3(...
    e(1,:), e(2,:), ...
    zeros(1, length(e)), 'g-', 'LineWidth', 2);
%Code for histagram
maxP = 0;
for side = 1:2
    % Calculate the histogram.
    p = [0 hist(r_xy(side,:), 20) 0];
    p = p / sum(p);
    maxP = max([maxP p]);
    dx = (axisMax(side) - axisMin(side)) / numel(p) / 2.3;
    p2 = [zeros(1,numel(p)); p; p; zeros(1,numel(p))]; p2 = p2(:);
    x = linspace(axisMin(side), axisMax(side), numel(p));
    x2 = [x-dx; x-dx; x+dx; x+dx]; x2 = max(min(x2(:), axisMax(side)), axisMin(side));

    % Calculate the curve.
    nPtsCurve = numel(p) * 10;
    xx = linspace(axisMin(side), axisMax(side), nPtsCurve);

    % Plot the curve and the histogram.
    if side == 1
        plot3(xx, ones(1, nPtsCurve) * axisMax(3 - side), spline(x,p,xx), 'r-', 'LineWidth', 2);
        plot3(x2, ones(numel(p2), 1) * axisMax(3 - side), p2, 'k-', 'LineWidth', 1);
    else
        plot3(ones(1, nPtsCurve) * axisMax(3 - side), xx, spline(x,p,xx), 'b-', 'LineWidth', 2);
        plot3(ones(numel(p2), 1) * axisMax(3 - side), x2, p2, 'k-', 'LineWidth', 1);
    end

end
grid on
hold off







