function [Cv_n] = calculate_Cv(dt,x,n)
v_initial = (x(2) - x(1))/dt;
v(1) = v_initial;
for i = 2:1:n
v(i) = (x(i) - x(i-1))/dt;
end
Cv_n = autocorr(v);
end


