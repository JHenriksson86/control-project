function [out] = make_init_traj(P0, P5)
%MAKE_INIT_TRAJ Summary of this function goes here
d_p = P5-P0;
out = [P0 (P0 + d_p*(1/5)) (P0 + d_p*(2/5)) (P0 + d_p*(3/5)) (P0 + d_p*(4/5)) P5]';

end

