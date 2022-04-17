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
[C,P] = blsprice(S(1:M),K, r, T, sigma);
subplot(2,1,1);
hold on;
plot(S(1:M), V(N, 1:M));
plot(S1, V1, '-o');
text(S1, V1, ['  ', 'Option price =', num2str(V1), ' at S = 100']);

title(['Option price using Implicit for M=', num2str(M), ', N=', num2str(N)]);
xlabel('Stock price'); ylabel('Option price at t = 0');

hold off;
subplot(2,1,2);
scatter(P, V(N,1:M), 10, 'filled');
refline;
title(['Compare blsprice vs Implicit method for ', num2str(M), ', N=', num2str(N)]);
xlabel('Option price using blsprice'); ylabel('Computed price');
