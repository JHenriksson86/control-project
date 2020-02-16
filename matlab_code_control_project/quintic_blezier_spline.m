
%% Intial trajectory

clear; close all;clc;

wp1 = [0 0]';
wp2 = [5 5]';
wp3 = [7 2]';
wp4 = [4 1]';

G1 = make_init_traj(wp1, wp2);
G2 = make_init_traj(wp2, wp3);
G3 = make_init_traj(wp3, wp4);

traj1 = zeros(2,10);
traj2 = zeros(2,10);
traj3 = zeros(2,10);
for i = 1:11
    t = (i-1)/10;
    traj1(:,i) = qb_spline(t, G1);
    traj2(:,i) = qb_spline(t, G2);
    traj3(:,i) = qb_spline(t, G3);
end

figure;
xlim([-1 8]);
ylim([-1 6]);
hold
sum_traj = [traj1 traj2 traj3];
plot(sum_traj(1,:), sum_traj(2,:), '--b');

%% Trajectory after first derivative heuristic

% Find tangent at wp2.
dt1 = derivative_tangent(wp1, wp2, wp3);

% Plot tangent
plot_tangent = [ wp2 (wp2 + dt1) ];
%plot(plot_tangent(1,:), plot_tangent(2,:), '-r');

%Find tangent at wp3
dt2 = derivative_tangent(wp2, wp3, wp4);

% Plot tangent
plot_tangent = [ wp3 (wp3 + dt2) ];
%plot(plot_tangent(1,:), plot_tangent(2,:), '-r');

% Set derivatives and calculate new trajectories
G1(5,:) = G1(6,:) - (1/5)*dt1';
G2(2,:) = (1/5)*dt1' + G2(1,:);

G2(5,:) = G2(6,:) - (1/5)*dt2';
G3(2,:) = (1/5)*dt2' + G3(1,:);

for i = 1:11
    t = (i-1)/10;
    traj1(:,i) = qb_spline(t, G1);
    traj2(:,i) = qb_spline(t, G2);
    traj3(:,i) = qb_spline(t, G3);
end
sum_traj = [traj1 traj2 traj3];
plot(sum_traj(1,:), sum_traj(2,:), '-.g');

%% Second derivatives heuristic

dtt2 = second_derivative_tangent(G1, G2);
G1(4,:) = (1/20)*dtt2 + 2*G1(5,:) - G1(6,:);
G2(3,:) = (1/20)*dtt2 + 2*G2(2,:) - G2(1,:);

dtt3 = second_derivative_tangent(G2, G3);
G2(4,:) = (1/20)*dtt3 + 2*G2(5,:) - G2(6,:);
G3(3,:) = (1/20)*dtt3 + 2*G3(2,:) - G3(1,:);

for i = 1:11
    t = (i-1)/10;
    traj1(:,i) = qb_spline(t, G1);
    traj2(:,i) = qb_spline(t, G2);
    traj3(:,i) = qb_spline(t, G3);
end
sum_traj = [traj1 traj2 traj3];
plot(sum_traj(1,:), sum_traj(2,:), '-r');

hold