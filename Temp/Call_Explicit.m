function [V0, V, alpha, beta] = Call_Explicit(S0, Smax, T, K, r,sigma, M,N)
%
%  Explicit method for solving the B-S equation of European call option
%
% Input
%   S0 -- current price of the stock
%   Smax -- right end of the stock price region [0, Smax]
%   T -- maturity of the option
%   K -- strike price
%   r -- risk free interest rate
%   simga -- volatility
%   M -- # of subintervals on [0, Smax]
%   N -- # of subintervals on [0, T]
%
% Output
%   V0 -- price for the option at t=0 and S0
%   V -- prices of the option at any t and S.
%   
%
%

%  determine dS and dtau to construct a grid of tau and S
dS = (Smax -0)/M;
dtau = T / N;
% initilize V, alpha and beta
V = zeros(N+1,M);
alpha= zeros(1,M);
beta = zeros(1,M);

% subinterval ends of S
S = 0: dS: Smax;

% initial condition of the PDE
for i = 1: M
    V(1, i) = max(S(i) - K, 0);    
end
% boundary condition on Smax
V(1, M) = Smax;
%
% main loop for the explicit method
%
for n = 1:N-1
    V(n+1, 1) = V(n, 1) * (1 - r * dtau); % boundary condition on S=0
    V(n+1, M) = V(n, M);   % boundary condition on Smax
    for i = 2:M-1
        %
        % Central differencing and up stream weighting
        %
        alpha_centeral = ((sigma ^ 2) * (S(i) ^2))/(((S(i) - S(i-1)) * (S(i+1) - S(i-1)))) - (r * S(i))/(S(i+1) - S(i-1));
        beta_centeral = ((sigma ^ 2) * (S(i) ^2))/(((S(i+1) - S(i)) * (S(i+1) - S(i-1)))) + (r * S(i))/(S(i+1) - S(i-1));
        if (alpha_centeral >= 0 ) && (beta_centeral >= 0)
            alpha(i) = alpha_centeral;
            beta(i) = beta_centeral;
        else 
            alpha_forward = ((sigma ^ 2) * (S(i) ^2))/(((S(i) - S(i-1)) * (S(i+1) - S(i-1))));
            beta_forward = ((sigma ^ 2) * (S(i) ^2))/(((S(i+1) - S(i)) * (S(i+1) - S(i-1)))) + (r * S(i))/(S(i+1) - S(i-1));
            alpha(i) = alpha_forward;
            beta(i) = beta_forward;
        end 
        %
        % updata V^{n+1} by V^n
        V(n+1, i) = V(n, i) * (1 - (alpha(i) + beta(i) + r) * dtau) + alpha(i) * dtau * V(n, i-1) + beta(i) * dtau * V(n, i+1);
    end  % end of M for-loop
end % end of N for-loop
%
% Evaluate the current option price V0
%
% find the smalleset interval including S0 
indx1 = max(find(S<=S0));
indx2 = min(find(S>=S0));
if indx1 == indx2  % S0 on the end of subintervals
    V0= V(N+1, indx1);
else    % S0 not on the end, estimate V0 by the linear interpolation
    w = (S0-S(indx1))/(S(indx2)-S(indx1));
    V0 = V(N+1,indx1)*w + (1-w)*V(N+1, indx2);
end      
end
    
