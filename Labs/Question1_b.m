clear all;
close all,
clc;
sigma = 0.3;
S1 = 20;
A1 = 0;
r = 0.05;
K = 20;
T = 0.25;
N = 30;
dT = T / N;
M = 100000;
P = zeros (1 , M);
for j = 1:M
    S = zeros (1, N);
    S(1) = 20;
    A = zeros (1, N);
    for i= 2:N
        p(i) = randn;
        S(i) = S(i-1) * exp ((r - (sigma.^2)/2)*dT + sigma * sqrt(dT) * p(i));
        A(i) = (1 / i)*sum (S);
    end 
    P(j) = max (A(N) - K, 0);
end 
V = exp (-r * T) * (1 / M) * sum (P)
