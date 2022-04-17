close all; clear all; clc;

% Basic input parameters
T = 1; % total time
S1 = 100; % Input current price
K = 100; % Strike price
r = 0.02; % Risk-free interest rate
sigma = 0.4; % Volatility

% The grid division factor (can take values 1,2,4,8 etc.)
% A higher c creates a finer grid
c = 4; 

% Create time grid using factor 'c'
dtau = (T/25)/c; % Time-step size

% Create price grid using factor 'c'
S = [0:0.1*K/c:0.4*K,...
0.425*K:0.05*K/c:0.8*K,...
0.805*K:0.02*K/c:0.9*K,...
0.905*K:0.01*K/c:1.1*K,...
1.12*K:0.02*K/c:1.2*K,...
1.25*K:.05*K/c:1.6*K,...
1.7*K:0.1*K/c:2*K,...
2.2*K, 2.3*K, 2.4*K, 2.6*K, 2.8*K,...
3.2*K, 3.4*K, 3.6*K, 4.3*K, 5*K, 7.5*K, 8.125*K, 8.75*K, 9.375*K, 10*K];

M = length(S); % Total steps along price axis
N = T/dtau; % Total steps along time axis

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
vector1 = [r*dtau, dtau *(alpha(1:M-1) + beta(1:M-1) + r), 0];
% vector2 and vector2 have an extra 0 each to support spdiags later below
vector2 = [0, 0, -dtau * beta(1:M-1)];
vector3 = [-dtau * alpha(1:M-1), 0, 0];
% Use spdiags to create a sparse matrix, to save memory
M_hat = spdiags(vector1', 0, M+1, M+1) ...
    + spdiags(vector2', 1, M+1, M+1) + spdiags(vector3', -1, M+1, M+1);
% Use LU decomposition for faster results
[L, U] = lu(I + M_hat);

for n = 1:N-1
    % Boundary conditions
    V(n+1, 1) = V(n, 1) / (1 - r * dtau);
    V(n+1, M) = V(n, M);
    % Solve using LU matrices calculated earlier for faster results
    B = transpose(V(n, :));
    V(n+1, :) = transpose(U\(L\B));
end % end of N for-loop

%
% Evaluate the current option price V1
%
% find the smalleset interval including S1 
indx1 = max(find(S<=S1));
indx2 = min(find(S>=S1));
if indx1 == indx2  % S1 on the end of subintervals
    V1= V(N, indx1);
else    % S1 not on the end, estimate V1 by the linear interpolation
    w = (S1-S(indx1))/(S(indx2)-S(indx1));
    V1 = V(N,indx1)*w + (1-w)*V(N, indx2);
end
disp(['Option price at (t=0) is ', num2str(V1), ' when S=', num2str(S1)]);

% Plot the graph and compare with blsprice 
hold on; [C,P] = blsprice(S(1:M),K, r, T, sigma); plot(S(1:M), P); plot(S(1:M), V(N, 1:M));
legend('Implicit method', 'blsprice');
xlabel('Stock price'); ylabel('Option price at t = 0');
