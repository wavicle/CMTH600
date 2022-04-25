close all;
clear all;
clc;
T = 0.25;
S1 = 0;
Smax = 3;
K = 1;
M = 100;
N = 20;
r = 0.05;
sigma = 0.3;
dS = (Smax)/M;
dtau = T / N;
V = zeros(N,M);
S = S1: dS: Smax;
for i = 1: M
    V(1, i) = max(S(i) - K, 0);    
end

for n = 1:N-1
    V(n+1, 1) = V(n, 1) * (1 - r * dtau);
    V(n+1, M) = V(n, M);
    for i = 2:M-1
        alpha_centeral(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i) - S(i-1)) * (S(i+1) - S(i-1)))) - (r * S(i))/(S(i+1) - S(i-1));
        beta_centeral(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i+1) - S(i)) * (S(i+1) - S(i-1)))) + (r * S(i))/(S(i+1) - S(i-1));
        alpha_forward(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i) - S(i-1)) * (S(i+1) - S(i-1))));
        beta_forward(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i+1) - S(i)) * (S(i+1) - S(i-1)))) + (r * S(i))/(S(i+1) - S(i-1));
        if (alpha_centeral(i) >= 0 ) && (beta_centeral(i) >= 0)
            
            alpha(i) = alpha_centeral(i);
            beta(i) = beta_centeral(i);
        else 

            alpha(i) = alpha_forward(i);
            beta(i) = beta_forward(i);
        end 
        V(n+1, i) = V(n, i) * (1 - (alpha(i) + beta(i) + r) * dtau) + alpha(i) * dtau * V(n, i-1) + beta(i) * dtau * V(n, i+1);
    end  % end of M for-loop
end % end of N for-loop
 V(N, M )
[Call, Put] = blsprice(Smax,K, r, T, sigma)
