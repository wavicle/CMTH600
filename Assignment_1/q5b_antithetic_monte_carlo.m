% Close old graphs, add random seed etc.
clear; clc; close all; hold off; rng(1);

% The various values of N for which we will run the simulation
N_values = [100, 500, 2000, 3000];

% The X and Y limits of the double integration
y_min = 0; y_max = 0.5;
x_min = 0; x_max = 2;

% The size of the domain for the limits of integration
m = (x_max - x_min) * (y_max - y_min);

% Stores the integral values based on simple Monte Carlo simulations
I_simple = zeros(1, length(N_values));
% Stores var(integral values) based on simple Monte Cancel simulations
I_simple_var = zeros(1, length(N_values));

% Stores the integral values based on antithetic variates method
I_anti = zeros(1, length(N_values));
% Stores the variance of integrals based on antithetic variates method
I_anti_var = zeros(1, length(N_values));

% Run simulation for various N values
for k = 1:length(N_values)
   N = N_values(k);
  
   % One vector for simple Monte Carlo, another for antithetic variates
   f_simple= zeros(1, N);
   f_anti= zeros(1, N);
   xsimple = zeros(1, N);
   ysimple = zeros(1, N);
   xanti = zeros(1, N);
   yanti = zeros(1, N);
  
   % We will average over N paths
   for i = 1:N
       % Compute y(i) and x(i) for the 'simple' Monte Carlo simulation
       ysimple(i)= y_min + (y_max - y_min)*rand;
       xsimple(i) = x_min + (x_max - x_min)*rand;
      
       % Compute y(i) and x(i) for the antithetic variates (see document
       % for detailed derivation)
       yanti(i) = y_min + y_max - ysimple(i);
       xanti(i) = x_min + x_max - xsimple(i);
      
       % Calculate for simple and antithetic variates
       f_simple(i) = xsimple(i)*ysimple(i)*ysimple(i);
       f_anti(i) = (f_simple(i) + (xanti(i)*yanti(i)*yanti(i)))/2;
   end
  
   % For both cases, calculate average over all paths
   I_simple(k) = m * sum(f_simple)/N;
   I_anti(k) = m * sum(f_anti)/N;
  
   % Also calculate the variances for comparison later
   I_simple_var(k) = var(f_simple);
   I_anti_var(k) = var(f_anti);
end

disp('The mean I values from simple Monte Carlo are:')
disp(I_simple);

disp('The var(I) values from simple Monte Carlo are:')
disp(I_simple_var);

disp('The mean I values from antithetic variates method are:')
disp(I_anti);

disp('The var(I) values from antithetic variates method are:')
disp(I_anti_var);

% Turn off scientific notation in Matlab (this is for plotting purposes)
ax = gca;
ax.YRuler.Exponent = 0;
hold on;

% Plot to compare variances in Simple vs Antithetic methods
plot(N_values, I_simple_var, '-O');
plot(N_values, I_anti_var, '-O');
legend('Simple Monte Carlo', 'Antithetic variates')
xlabel('N');
ylabel('var(integral)');
title('Simple and Antithetic variances of the integral');
