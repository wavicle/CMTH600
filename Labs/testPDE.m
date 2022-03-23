clc; clear;
S0=100; 
Smax = 3;
T = 0.25;
K= 1;
r =0.05;
sigma=0.3;
M = 200;
N = 200;
% explicit method
%[V0, V] = Call_explicit(S0, Smax, T, K, r,sigma, M,N);
% implicit method
VV = PDE_Implicit(S0, Smax, T, K, r,sigma, M,N);
% B-S formula
[C,P] = blsprice(S0, K, r, T, sigma);
