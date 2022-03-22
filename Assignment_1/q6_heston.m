clear; clc; close all; rng(1); % Set seed for repeatable samples
hold on;

M = 10; % Number of simulations
r = 0.02; % Interest rate
dT = 1/100; % The timestep size delta T
T = 1; % expiry time
N = T/dT; % Number of timesteps
df = exp(-r*T); % Discount factor considering continuous compounding

v = zeros(1, N); % Variance at each timestep
v(1) = 0.0174; % Initial variance (Matlab index starts with 1, not 0)
v_bar = 0.0354; % Long term mean of variance

eta = 0.3877; % The value of η in the problem
rho = -0.7165; % The correlation coefficient ρ
lambda = 1.3253; % The lambda value

K = 100; % Strike price

x = zeros(1, N); % Ln(stock price) for each time step
S = zeros(1, N); % Stock price at each time step
S(1) = K; % Matlab index starts with 1 and this is an at-the-money call
x(1) = log(S(1)); % Initial value of x(i)

% Store value for each path in f which we will average later
f = zeros(1, M);

% Run simulation over M paths
for path = 1: M
  
   % phi1 and phi2 are ~ N(0,1) with correlation = rho
   phi1 = zeros(1, N);
   phi2 = zeros(1, N);

   % Calculate S at N timesteps. Since we are calculating the (i + 1)th
   % value, the loop runs from 1 to N - 1
   for i = 1: N - 1
       % phi1 and phi2 are ~ N(0,1) with correlation = rho
       phi1(i) = randn;
       phi2(i) = (rho*phi1(i)) + (sqrt(1 - rho*rho)*randn);
      
       x(i+1) = x(i) + (r*dT) - (v(i)*dT/2) + sqrt(v(i)*dT)*phi1(i);
       S(i+1) = exp(x(i));
      
       % Long equation, let's divide into three parts for simplicity
       part1 = ( sqrt(v(i)) + (eta/2)*sqrt(dT)*phi2(i) )^2;
       part2 = -lambda*(v(i) - v_bar)*dT;
       part3 = -(eta*eta/4)*dT;
      
       % Make sure v(i) is non-negative
       v(i + 1) = abs(part1 +part2 + part3);
   end
  
   % Value = Discounted payoff for European call
   f(path) = df*max(S(N) - K, 0);
   plot(S);
end

xlabel('Time'); ylabel('Stock price'); title('Simulation of 10 sample paths under Heston''s model');

% Option price is expected discounted payoff
V0 = sum(f)/M;
disp('The option price is: '); disp(V0);
