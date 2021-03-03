function [x_without_inertia] = distance_no_inertia(dt, wi, D, n)

x_without_inertia(1) = 0;
    for i = 2:n
    x_without_inertia(i)=x_without_inertia(i-1)+(sqrt(2 * D * dt)) *wi(i);
    end 
end