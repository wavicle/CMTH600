% Use rng to fix random numbers
clc; clear; close all; rng(54);
disp('Please give this program a few min to complete')

% Use values from previous questions
S0 = 1;
T = 0.25; 
r = 0.03; 
M = 10000; 
N = 100;

% Given market values
K = [0.80; 0.85; 0.90; 0.95; 1.00; 1.05; 1.10];
V0_mkt = [0.3570; 0.2792; 0.2146; 0.1747; 0.1425; 0.1206; 0.0676];

% Starting point vectors
x0_1  = [0.2; 0.0; 0.0];
x0_2 = [0.2; 0.1; 0.01];

options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true);
options.Algorithm = 'levenberg-marquardt';
options.Display = 'iter';

% Calibrate with both starting vectors and display the outputs
[x1, err_1] = lsqnonlin(@myfun, x0_1, [], [], options);
[x2, err_2] = lsqnonlin(@myfun, x0_1, [], [], options);
disp(['Best Fit X1 = [', num2str(x1'), '] with error = ', num2str(err_1)]);
disp(['Best Fit X2 = [', num2str(x2'), '] with error = ', num2str(err_2)]);

% Calculate the option values with best fits and compare with market
calc_x1 = @(k_val) Eur_Call_LVF_MC(S0, k_val, T, r, x1, M, N);
calc_x2 = @(k_val) Eur_Call_LVF_MC(S0, k_val, T, r, x2, M, N);
V0_x1 = arrayfun(calc_x1, K); V0_x2 = arrayfun(calc_x2, K);
hold off; figure(1); hold on;
plot(K, V0_mkt, '-o'); plot(K, V0_x1, '-o'); plot(K, V0_x2, '-o');
legend('Market Price', 'MC with best fit 1', 'MC with best fit 2');
title('Option price comparison');
xlabel('Strike Price K'); ylabel('Option Price');

% Calculate and compare implied volatities
impvol_mkt = blsimpv(S0, K, r, T, V0_mkt);
impvol_x1 = blsimpv(S0, K, r, T, V0_x1);
impvol_x2 = blsimpv(S0, K, r, T, V0_x2);
hold off; figure(2); hold on;
plot(K, impvol_mkt, '-o'); plot(K, impvol_x1, '-o'); plot(K, impvol_x2, '-o');
legend('From Market Price', 'From best fit 1', 'From best fit 2');
title('Implied volatility comparison');
xlabel('Strike Price K'); ylabel('Implied Volatility');


