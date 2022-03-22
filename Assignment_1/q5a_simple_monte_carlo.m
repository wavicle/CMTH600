% Clear screen, set random seed etc.
clear; clc; close all; rng(100);

% The various values of N for which we will run the simulation
N_values = [100, 500, 2000, 3000];

% The X and Y limits of the double integration
y_min = 0; y_max = 0.5;
x_min = 0; x_max = 2;

% The size of the domain for the limits of integration
m = (x_max - x_min) * (y_max - y_min);

% Stores the integral values based on simple Monte Carlo simulations
I_simple = zeros(1, length(N_values));

% Stores var(integral values) based on simple Monte Carlo simulations
I_simple_var = zeros(1, length(N_values));

% Run simulation for various N values
for k = 1:length(N_values)
   N = N_values(k);
  
   % The values calculated for each path
   f_simple= zeros(1, N);
   x = zeros(1, N);
   y = zeros(1, N);
  
   % We will average over N paths
   for i = 1:N
       % Compute y(i) and x(i) using shifted and scaled U(0,1)
       y(i) = y_min + (y_max - y_min)*rand;
       x(i) = x_min + (x_max - x_min)*rand;
      
       % Calculate f(i) for the simple Monte Carlo method
       f_simple(i) = x(i)* y(i)* y(i);
   end
  
   % Calculate average over all paths
   I_simple(k) = m * sum(f_simple)/N;
   % Also calculate the variances for comparison later
   I_simple_var(k) = var(f_simple);
end

% Print the various values on the console
disp('Mean I values for N=100,500,2000,3000 from simple Monte Carlo are:')
disp(I_simple);

disp('var(I) N=100,500,2000,3000 from simple Monte Carlo are:')
disp(I_simple_var);
