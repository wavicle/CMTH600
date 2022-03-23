close all; clear all; clc;

% Basic input parameters
T = 1; % total time
S1 = 100; % Input current price
K = 100; % Strike price
r = 0.02; % Risk-free interest rate
sigma = 0.4; % Volatility
dtau = 0.005; % Time-step size

% Choose smax using 3-sigma rule as explained in the class
Smax = S1*exp((r-0.5*sigma*sigma)*T + 3*sigma*sqrt(T));

% Create price grid
S = [0:0.1*K:0.4*K,...
0.45*K:0.05*K:0.8*K,...
0.82*K:0.02*K:0.9*K,...
0.91*K:0.01*K:1.1*K,...
1.12*K:0.02*K:1.2*K,...
1.25*K:.05*K:1.6*K,...
1.7*K:0.1*K:2*K,...
2.2*K, 2.4*K, 2.8*K,...
3.6*K, 5*K, 7.5*K, 10*K];

M = length(S); % Total steps along price axis
N = T/dtau; % Total steps along time axis

% Initialize option price V(tau, S) with initial and boundary conditions
V = zeros(N,M);
for i = 1: M
    V(1, i) = max(K - S(i), 0); % Payoff for PUT option at tau = 0
end
V(1 , M) = 0; % Boundary condition at Smax for PUT option

% Initialize some matrices for later use
alpha = zeros(1,M);
beta = zeros(1,M);
I = eye(M);

for n = 1:N-1
    V(n+1, 1) = V(n, 1) * ((1 + (r / 2) * dtau) / (1 - (r / 2) * dtau));
    V(n+1, M) = V(n, M);
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
    vector1 = [-r*dtau/2, dtau/2 *(alpha(1:M-2) + beta(1:M-2) + r), 0];
    vector2 = [0, -dtau/2 * beta(1:M-2)];
    vector3 = [-dtau/2 * alpha(1:M-2), 0];
    M_hat = diag(vector1) + diag(vector2, 1) + diag(vector3, -1);

    V(n+1, :) = transpose(inv(I + M_hat) *transpose(V(n, :))) ...
        * transpose(I - M_hat);
end % end of N for-loop
