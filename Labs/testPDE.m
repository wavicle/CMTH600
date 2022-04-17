clc; clear;
S0=100; 
Smax = 10000;
T = 1;
K= 100;
r =0.02;
sigma=0.4;
M = 400;
N = 400;
% explicit method
[V0, V] = Call_explicit(S0, Smax, T, K, r,sigma, M,N);
% implicit method
[V0,VV] = PDE_Implicit(S0, Smax, T, K, r,sigma, M,N);
% B-S formula
[C,P] = blsprice(S0, K, r, T, sigma);
P