function [MSD] = mean_square_displacement(x,n)
numberOfDeltat = floor(sqrt(n));
MSD = zeros(numberOfDeltat,1);
for i = 1:numberOfDeltat
  MSD(i,1) =mean(sum((x(1+i:end)-x(1:end-i)).^2));
end

end