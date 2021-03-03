function [r_x,r_y] = tweezer_position_XYplane(dt,k_xy, wi, D, g, n)
 r_xy = ones(2,n);
 r_xy(:,1) = [0;0];
 for i = 2:n
    r_xy(:,i) = r_xy(:,i-1) - (1/g)*k_xy*r_xy(:,i-1)*dt + sqrt(2*D*dt) * wi(:,i);
    r_x(i) = r_xy(1,i);
    r_y(i) = r_xy(2,i); %separating vector r_xy into r_x and r_y
 end

end