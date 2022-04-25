function V = Eur_Call_LVF_MC(S0, K, T, r, x, M, N)
%
% Price the European call option of the LVF model by the Monte Carlo method
%
% Input:
% S0 - initial stock price
% K  - strike price
% T  - time to maturity
% r  - risk free interest rate
% x  - column vector parameters for the LVF sigma, [x1; x2; x3]
% M  - number of simulated paths
% N  - number of time steps, i.e., dt = T/N
%
% Output:
% V - European call option price at t = 0 and S0.

dT = T/N; % The timestep size delta T
logS = zeros(1, N); % Ln(stock price) for each time step
S = zeros(1, N); % Stock price at each time step
S(1) = S0; % Start with S0 for the price
logS(1) = log(S(1)); % Initial value of x(i)

% Store value for each path in f which we will average later
f = zeros(1, M);

% Run simulation over M paths
for path = 1: M
   % phi is a normal variable: phi ~ N(0,1)
   phi = zeros(1, N);

   % Calculate S at N timesteps. Since we are calculating the (i + 1)th
   % value, the loop runs from 1 to N - 1
   for i = 1: N - 1
       phi(i) = randn;
       % sigma is time-dependent and needs to be calculated at each step
       sigma = max(0, [1 S(i) S(i)*S(i)]*x);
       logS(i+1) = logS(i) + (r*dT) - (sigma*sigma*dT/2) + sigma*sqrt(dT)*phi(i);
       S(i+1) = exp(logS(i));
   end
  
   % Value = Payoff for European call
   f(path) = max(S(N) - K, 0);
end

% Option price is average discounted payoff
V = exp(-r*T)*sum(f)/M;
end