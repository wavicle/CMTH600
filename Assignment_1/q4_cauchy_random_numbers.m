% Summary: Since we are not given the value of c, we will use c = 7 as an
% example. We will generate 1000 numbers and also compare the generated
% CDF to the actual CDF (from formula).
clear; clc; close all; hold on; rng(10); % Clear everything

c = 7; % c = 7 is just an example since c is not specified in assignment
total = 1000; % We will generate 10000 random numbers
results = zeros(1, total); % This vector stores generated random values
U = zeros(1, total);

% Generate Cauchy numbers in a loop
for i = 1:total
    % Generate a random number ~ U(0,1)
    U(i) = rand;
    results(i) = c * tan(pi* (U(i) - 0.5)); % F_inverse (U) as derived
end

% Print the results
disp("Mean of random numbers: "); disp(mean(results));
disp("Variance of random numbers: "); disp(var(results));
% Plot the CDF of this distribution
cdfplot(results);

% Next we superimpose the actual CDF on this for comparison
x = min(results):0.1:max(results); % Same range as empirical CDF
plot(x, (1/pi)*atan(x/c) + 0.5);
legend('Generated Cauchy CDF', 'Actual Cauchy CDF')
title('Graph comparing generated vs actual Cauchy CDF for c = 7');