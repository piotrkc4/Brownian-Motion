function [Var] = averaging(x1,mean)
Va = ones(length(x1),1)
Va(1) = (x1(1) - mean)^2;
Var(1) = Va(1);
for i = 2:length(x1) %averaging over x
    Va(i) = (x1(i) - mean)^2;
    Var(i) = (Va(i) + Va(i-1))/i
end