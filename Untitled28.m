clc
clear all
close all
w = zeros(100,10000);
M=zeros(100,1);
V=zeros(100,1);
for i=1:100
    w(i,:) = randn(10000,1);
M(i) = mean(w(i,:));
V(i) = var(w(i,:));
iter(i) = i;
end
figure
scatter(iter,M)
ylabel('mean')
xlabel('iteration')
title('Means of distributions generated by randn')
figure
scatter(iter,V)
ylabel('variance')
xlabel('iteration')
title('Variances of distributions generated by randn')