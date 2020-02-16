function [tangent] = derivative_tangent(old_wp_sum,wp1, wp2)
%DERIVATIVE_TANGENT Summary of this function goes here

w1 = wp1 - old_wp_sum;
w2 = wp2 - old_wp_sum;
magnitude = min(norm(w1), norm(w2-w1))/2;
tangent = magnitude * (w1 + (w2-w1))/norm(w1 + (w2-w1));

end

