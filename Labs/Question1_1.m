clear all;
close all,
clc;
sigma = 0.3;
S0 = 20;
A1 = 0;
r = 0.05;
K = 20;
T = 0.25;
N = [1000 10000 100000];
for j= 1:length(N)
    S_T = zeros(1 , N(j));
    f = zeros(1 , N(j));
    for i = 1:N(j)
        S_T(i) = S0 * exp ((r - (sigma.^2)/2)* T + sigma * sqrt(T) * randn);
        f(i) = max(S_T(i) - K, 0);
    end
    V = exp (-r * T) * (1 / N(j)) * sum (f)
end     