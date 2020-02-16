function [tangent] = second_derivative_tangent(spline1, spline2)
%SECOND_DERIVATIVE_TANGENT Summary of this function goes here
wp1 = spline1(1,:);
wp2 = spline2(1,:);
wp3 = spline2(6,:);

t1 = 5*(spline1(2,:) - spline1(1,:));
t2 = 5*(spline2(2,:) - spline2(1,:));
t3 = 5*(spline2(6,:) - spline2(5,:));

d12 = norm(wp2-wp1);
d23 = norm(wp3-(wp2-wp1));

alpha = d23 / (d12 + d23);
beta = d12 / (d12 + d23);
 
tangent = alpha*(6*wp1+2*t1+4*t2-6*wp2)+beta*(-6*wp2-4*t2-2*t3+6*wp3);

end

