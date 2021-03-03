figure
plot3(xt1, yt1, zeros(size(r_xyz(1:2,:))), 'k.', 'MarkerSize', 1);
xlabel('xt (nm)')
ylabel('yt (nm)')
M = mean(r_xyz(1:2,:));
data_zeroMean = bsxfun(@minus, r_xyz(1:2,:), M);
[V,D] = eig(data_zeroMean * data_zeroMean' / (size(data_zeroMean, 1))');
[D, order] = sort(diag(D), 'descend');
D = diag(D);
V = V(:, order);
V = V * sqrt(D);
t = linspace(0, 2 * pi);
e = bsxfun(@plus, 2*V * [cos(t); sin(t)], M');
plot3(...
    e(1,:), e(2,:), ...
    zeros(1, length(e)), 'g-', 'LineWidth', 2);