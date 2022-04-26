function V0 = Eur_Call_LVF_PDE(S0, K, T, r, x, SMax, M, N)
%
% Price the European call option of the LVF model by the explicit Finite difference method.
%
% Input
% S0 	- initial stock price
% K 	- strike price
% T 	- maturity
% r 	- risk free interest rate
% x 	- vector parameters for the LVF sigma, [x1; x2; x3]
% Smax 	- upper bound of the stock price
% M 	- number of stock price difference, i.e., dS = Smax/M
% N 	- number of time steps, i.e., dT = dTau = T/N.
%
% Output
% V0 	- European call option price at t = 0 and S0.

% Initialize the price vector
dS = SMax/M;
S = 0:dS:SMax;
dTau = T/N;

% Initialize some matrices for later use
alpha = zeros(1,M);
beta = zeros(1,M);

% Initialize option price V(tau, S) with initial and boundary conditions
V = zeros(N, M);
for i = 1: M
    V(1, i) = max(S(i) - K, 0); % Payoff for CALL option at tau = 0
end

% Iterate over the asset price
for n = 1:N-1
    % Boundary conditions
    V(n+1, 1) = V(n, 1) * (1 - r * dTau);
    V(n+1, M) = V(n, M);
    
    % Iterate over the time to maturity
    for i = 2:M-1
        % sigma is price-dependent and needs to be calculated at each step
        sigma = max(0, [1 S(i) S(i)*S(i)]*x);
        
        % Find alpha, beta with central difference, upstream weighting
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
        
        % Calculate the next value of V from previous values
        V(n + 1, i) = V(n , i)*(1 - (alpha(i) + beta(i) + r)*dTau) ...
            + alpha(i)*dTau*V(n, i - 1)...
            + beta(i)*dTau*V(n, i + 1);
    end
end % end of N for-loop

%
% Evaluate the current option price V0 using interpolation
%
% find the smalleset interval including S0 
indx1 = max(find(S<=S0));
indx2 = min(find(S>=S0));
if indx1 == indx2  % S0 on the end of subintervals
    V0= V(N, indx1);
else    % S0 not on the end, estimate V1 by the linear interpolation
    w = (S0-S(indx1))/(S(indx2)-S(indx1));
    V0 = V(N,indx1)*w + (1-w)*V(N, indx2);
end

end


