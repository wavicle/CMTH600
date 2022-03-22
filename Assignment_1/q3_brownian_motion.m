% Clear all existing Matlab variables
clear; clc; close all; rng(2000);

% This is the number of paths used in the simulation
total_paths = 100000;

% The time horizon (T) for which we are running the simulation
T = 1;

% Values of various timesteps used in this simulation
K = [100, 200, 400, 800, 1600];

% Use vectorized way to get Δt = T/K (we will use it for plotting graphs)
dT_values = T./K;

% Used to store the "mean" value for each value of K
I_mean = zeros(1, length(K));

% Used to store var(I) for each value of K
I_var = zeros(1, length(K));

% Summary: We have to simulate for various values of K and
% each simulation has to run for 100000 paths. So the outer for loop
% iterates over the K values and the inner for loop runs over the paths.
% We now simulate I for various values of K (total timesteps)
for i = 1: length(K)
  
   % Δt = current time-step-size
   dT = T/K(i);
  
   % Initialize the path valus, which will store the I calculated for 
   % each path (used for variance comparison later). 
   path_values = zeros(1, total_paths);   
  
   % Law of Large Numbers: To get I(Δt) we simulate over 100000 paths
   % and then take the average at the end of the loop
   for path = 1: total_paths
       % Use vectorization instead of for-loops for smaller code. This is
       % the same as generating N(0,1) from 1 to K.
       phi = normrnd(0, 1, [1, K(i)]);
       % For each path use given formula: path_value(path) = ∑ Φ^2(i)*Δt
       path_values(path) = sum((phi.^2)*dT);
   end
  
   % The I(Δt) is the average along all paths for current K
   I_mean (i) = sum(path_values)/total_paths;
   % We store the variance of I(Δt) for current K for plotting later
   I_var (i) = var(path_values);
end

% Print I(Δt) and var(I) for each Δt on the console
for i = 1: length(K)
   fprintf('When Δt = %f, I(Δt) = %f and var(I) = %f \n', dT_values(i), I_mean(i), I_var(i));
end

% Let's plot var(I) vs Δt - we expect this to be a straight line
figure
plot(dT_values, I_var, '-O');
title('How var(I) changes with Δt')
xlabel('Values of Δt')
ylabel('Variance of I(Δt)')
