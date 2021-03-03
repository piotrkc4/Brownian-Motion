function [Cx] = calculate_Cv_in3D(dt,x,n)
v_initial = (x(:,2) - x(:,1))/dt;
v(:,1) = v_initial;
for i = 2:1:n
v(:,i) = (r_xyz(:,i) - r_xyz(:,i-1))/dt;
end
Cx = autocorr(v);
end