clc
clear all
close all
for i = 1:100
wi(:,i) = randn(100,1);
end
for i = 1:100
Average(i) = mean(wi(:,i));
Variance(i) = var(wi(:,i));
num_iter(i) = i 
end
figure
scatter(num_iter,Average)
xlim([0 100])
figure
scatter(num_iter,Variance,10)