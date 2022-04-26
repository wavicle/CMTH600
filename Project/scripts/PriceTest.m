clc; clear; close all;

% Result from Monte-Carlo simulation
V_MC = Eur_Call_LVF_MC(1, 1, 0.25, 0.03, [0.2; 0.003; 0.001], 10000, 100);
disp(['Monte-Carlo result = ', num2str(V_MC)]);

% Result from Explicit Finite Difference method
V_PDE = Eur_Call_LVF_PDE(1, 1, 0.25, 0.03, [0.2; 0.003; 0.001], 3, 30, 100);
disp(['Explicit Finite-Diff result = ', num2str(V_PDE)]);
