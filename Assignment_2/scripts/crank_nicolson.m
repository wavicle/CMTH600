close all; clear all; clc;

% Basic input parameters
T = 1; % total time
S1 = 100; % Input current price
K = 100; % Strike price
r = 0.02; % Risk-free interest rate
sigma = 0.4; % Volatility
M = 800; % Total steps along price axis
N = 800; % Total steps along time axis

% Find Smax using the 3-sigma rule, then create the price grid
S0 = 0;
SMax = S1*exp((r-0.5*sigma*sigma)*T + 3*sigma*sqrt(T));
dS = (SMax - S0)/M;
S = S0:dS:SMax;
dtau = T/N;

% Initialize option price V(tau, S) with initial and boundary conditions
V = zeros(N,M + 1);
for i = 1: M
    V(1, i) = max(K - S(i), 0); % Payoff for PUT option at tau = 0
end
V(1 , M+1) = 0; % Boundary condition at Smax for PUT option

% Initialize some matrices for later use
alpha = zeros(1,M+1);
beta = zeros(1,M+1);
I = eye(M+1);

% Calculate M_hat (only once) and the LU decomposition of (I + M_hat). 
% As per the assignment instructions, we should avoid unnecessary
% calculations, hence this is pulled outside the second loop.
for i = 2:M-1
    alpha_central(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i) - S(i-1)) ...
        * (S(i+1) - S(i-1)))) - (r * S(i))/(S(i+1) - S(i-1));
    beta_central(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i+1) - S(i)) ...
        * (S(i+1) - S(i-1)))) + (r * S(i))/(S(i+1) - S(i-1));
    alpha_forward(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i) - S(i-1)) ...
        * (S(i+1) - S(i-1))));
    beta_forward(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i+1) - S(i)) ...
        * (S(i+1) - S(i-1)))) + (r * S(i))/(S(i+1) - S(i-1));
    if (alpha_central(i) >= 0 ) && (beta_central(i) >= 0)
        alpha(i) = alpha_central(i);
        beta(i) = beta_central(i);
    else 
        alpha(i) = alpha_forward(i);
        beta(i) = beta_forward(i);
    end 
end  % end of M for-loop
vector1 = [-r*dtau/2, dtau/2 *(alpha(1:M-1) + beta(1:M-1) + r), 0];
% vector2 and vector2 have an extra 0 each to support spdiags later below
vector2 = [0, 0, -dtau/2 * beta(1:M-1)];
vector3 = [-dtau/2 * alpha(1:M-1), 0, 0];
% Use spdiags to create a sparse matrix, to save memory
M_hat = spdiags(vector1', 0, M+1, M+1) ...
    + spdiags(vector2', 1, M+1, M+1) + spdiags(vector3', -1, M+1, M+1);
% Use LU decomposition for faster results
[L, U] = lu(I + M_hat);

for n = 1:N-1
    % Boundary conditions
    V(n+1, 1) = V(n, 1) * ((1 + (r / 2) * dtau) / (1 - (r / 2) * dtau));
    V(n+1, M) = V(n, M);
    % Solve using LU matrices calculated earlier for faster results
    B = (I - M_hat)*transpose(V(n, :));
    V(n+1, :) = transpose(U\(L\B));
end % end of N for-loop
