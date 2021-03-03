function [r_xyz] = tweezer_position(dt,k, wi, D, g, n)
r_xyz(:,1) = [0;0;0];
for i = 2:n
    r_xyz(:,i) = r_xyz(:,i-1) - (1/g)*k*r_xyz(:,i-1)*dt + sqrt(2*D*dt) * wi(:,i);
end
r_xyz = r_xyz/1e-09; %convert to nanometers
end