close all;
clear all;
clc;
T = 0.25;
S1 = 0;
Smax = 3;
K = 1;
M = 200;
N = 200;
r = 0.05;
sigma = 0.3;
dS = (Smax - 0)/M;
dtau = T / N;
V = zeros(N,M);
V(1 , M) = Smax;
S = 0: dS: Smax;
for i = 1: M
    V(1, i) = max(S(i) - K, 0);    
end
alpha = zeros(1,M);
beta = zeros(1,M);
for n = 1:N-1
    V(n+1, 1) = V(n, 1) * ((1 + (r / 2) * dtau) / (1 - (r / 2) * dtau));
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
           end  % end of M for-loop
    vector1 = [dtau/2 *(alpha(1:M-1) + beta(1:M-1) - r), 0];
    vector2 = [0, -dtau/2 * beta(1:M-2)];
    vector3 = [-dtau/2 * alpha(1:M-2), 0];
    M_hat = diag(vector1) + diag(vector2, 1) + diag(vector3, -1);
    I = eye(M);
     V(n+1, :) = transpose(inv(I + M_hat) *transpose(V(n, :))) * transpose(I - M_hat);
end % end of N for-loop
V(N, M-1)
%[Call, Put] = blsprice(Smax,K, r, T, sigma)

            
    
