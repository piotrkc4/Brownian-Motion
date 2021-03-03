function [x_inertia] = distance_with_inertia(dt, wi, g, m, T, kb, n)
F = g/m;
M = dt * F;
A = sqrt(2 * kb * T * g);
B = (dt)^(3/2);
C = 1 + M;
D = 2 + M;
E = m * C;
%x_inertia = ones(n,1);
x_inertia(1) = 0;
x_inertia(2) = (A/E) * B * wi(2);
    for i = 3:n
        Term_1 = (D/C)*x_inertia(i-1);
        Term_2 = (1/C) * x_inertia(i-2);
        Term_3 = (A/E) * B * wi(i);
        x_inertia(i) = Term_1 - Term_2 + Term_3;
    end 
end