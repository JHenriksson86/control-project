clear; close all; clc

% Control variables
syms v real;
syms w real;

% State variables
syms x real;
syms y real;
syms theta real;

syms x_next real;
syms y_next real;
syms theta_next real;

syms x_dot real;
syms y_dot real;
syms theta_dot real;

% Constants
syms dt real;
syms r real;
syms L real;

A = [0 0 -r*v*sin(theta); 0 0 r*v*cos(theta); 0 0 0];
X = [x y theta]';

B = [r*cos(theta) 0; r*sin(theta) 0; 0 r/L];
u = [v w]';

%x_dot = (x_next-x)/dt;
%y_dot = (y_next-y)/dt;
%theta_dot = (theta_next-theta)/dt;
X_dot = [x_dot; y_dot; theta_dot];

u_sum = pinv(B)*(X_dot - A*X)