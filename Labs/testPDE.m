S0=100; 
Smax = 500;
T = 0.5;
K= 100;
r =0.02;
sigma=0.4;
M = 200;
N=1000;
dtau=T/N;
dS = Smax/M;
dtau<(dS)^2/(sigma^2*Smax^2)
% explicit method
[V0, V] = Call_explicit(S0, Smax, T, K, r,sigma, M,N);
% implicit method
[V1, VV] = PDE_Implicit(S0, Smax, T, K, r,sigma, M,N);
% B-S formula
[C,P] = blsprice(S0, K, r, T, sigma)