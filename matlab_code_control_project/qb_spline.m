function [result] = qb_spline(t, G)
%QB_SPLINE Quintic Bezier Spline function
% t should be between 1 and 0
% G = [P0 P1 P2 P3 P4 P5 P6]'

T = [t^5 t^4 t^3 t^2 t 1];
M = [-1 5 -10 10 -5 1;
     5 -20 30 -20 5 0;
     -10 30 -30 10 0 0;
     10 -20 10 0 0 0;
     -5 5 0 0 0 0;
     1 0 0 0 0 0];

result = T * M * G;
end

