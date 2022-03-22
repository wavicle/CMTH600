clear all;
close all,
clc;
sigma = 0.3;
S1 = 20;
A1 = 0;
r = 0.05;
K = 20;
dT = [0.1, 0.01, 0.005, 0.001];
T = 1;
for j=1:length(dT)
    N = (1/dT(j));
    M = (1/dT(j))^2;
    P = zeros (1 , M);
    for i = 1:M
       S = zeros (1, N);
       A = zeros (1, N);
       S(1) = S1;
       A(1) = S(1);
       for q= 2:N
          p = randn;
          S(q) = S(q-1) * exp ((r - (sigma.^2)/2)*dT(j) + sigma * sqrt(dT(j)) * p );
          A(q) = (1 / q) * sum(S);
      end 
      P(i) = max (A(N) - K, 0);
    end 
    V(j) = exp (-r * T) * (1 / M) * sum (P)
end
semilogx(dT, V, '-*'); grid;
% plot(dT, V, '-*'); grid;
set(gca, 'XDir','reverse')
xlabel('\DeltaT'),
ylabel('V')
