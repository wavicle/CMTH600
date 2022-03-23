close all; clear; clc;
% This can be 'CrankNicolson' or 'Implicit'
method = 'CrankNicolson';
% Option type can be PUT or CALL
optionType = 'PUT';
% Input params as given in the assignment
T = 1;
S1 = 100;
K = 100;
M = 500;
N = 100;
r = 0.02;
sigma = 0.4;
% Choose smax using 3-sigma rule as explained in the class
Smax = S1*exp((r-0.5*sigma*sigma)*T + 3*sigma*sqrt(T));
dS = (Smax)/M;
dtau = T/N;
I = eye(M);
V = zeros(N, M);
S = 0: dS: Smax;
for i = 1:M
    switch optionType
    case 'PUT'
        V(1, i) = max(K - S(i), 0);
    case 'CALL'
        V(1, i) = max(S(i) - K, 0);
    end
end
alpha = zeros(1, M);
beta = zeros(1, M);
% Iterate over all time upto T
for n = 1:N-1
% We have combined both methods using 'theta' as taught in Week 8
switch method
    case 'CrankNicolson'
        V(n+1, 1) = V(n, 1) * ((1 + (r/2)*dtau)/(1-(r/2)*dtau));
        theta = 0.5;
    case 'Implicit'
        V(n+1, 1) = V(n, 1) / (1 + r * dtau);
        theta = 0;
end
% Iterate over all possible S
for i = 2:M-1
    alpha_centeral(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i) - S(i-1)) ...
    * (S(i+1) - S(i-1)))) - (r * S(i))/(S(i+1) - S(i-1));
    beta_centeral(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i+1) - S(i)) ...
    * (S(i+1) - S(i-1)))) + (r * S(i))/(S(i+1) - S(i-1));
    alpha_forward(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i) - S(i-1)) ...
    * (S(i+1) - S(i-1))));
    beta_forward(i) = ((sigma ^ 2) * (S(i) ^2))/(((S(i+1) - S(i)) ...
    * (S(i+1) - S(i-1)))) + (r * S(i))/(S(i+1) - S(i-1));
    % Give preference to central params, else choose forward param
    % since it will always be positive (upstream weighting)
    if (alpha_centeral(i) >= 0 ) && (beta_centeral(i) >= 0)
        alpha(i) = alpha_centeral(i);
        beta(i) = beta_centeral(i);
    else
        alpha(i) = alpha_forward(i);
        beta(i) = beta_forward(i);
    end
end
% We use spdiags instead of full matrix objects for efficiency
vectorMid = [dtau/2*(alpha(1:M-1) + beta(1:M-1) - r), 0];
vectorUp = [0, 0, -dtau/2*beta(1:M-2)];
vectorLow = [-dtau/2*alpha(1:M-2), 0, 0];
M_hat = spdiags([vectorLow' vectorMid' vectorUp'], -1:1, M, M);
% Perform an LU decomposition to solve for V(n+1)
MCommon = 2*M_hat;
[Lower, Upper] = lu(I + (1 - theta)*MCommon);
B = (I - theta*MCommon)*V(n, :)';
Y = Lower\B;

V(n+1, :) = Upper\Y;
end
% Plot the values obtained at tau = T (t = 0)
plot(S(1:M), V(N, :));
xlabel('Stock price'); ylabel('Option value');