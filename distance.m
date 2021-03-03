function [x] = distance(n, wi, dt)
x = zeros(n,1); %define x as vertical array of size (n,1)
    for i = 2:n
    x(i)=x(i-1) + sqrt(dt)*wi(i); %from paper
    end
end