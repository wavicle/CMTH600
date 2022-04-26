function [F, J] = myfun(x)
% Use values from previous questions
S0 = 1;
T = 0.25; 
r = 0.03; 
M = 10000; 
N = 100;

% The differential change in x
dx = [0.01, 0.01, 0.01];

% Strike prices from market information
K = [0.80; 0.85; 0.90; 0.95; 1.00; 1.05; 1.10];
V0_mkt = [0.3570; 0.2792; 0.2146; 0.1747; 0.1425; 0.1206; 0.0676];

% Define this as an inner function for repeated use
function val = Fx(x)
    val = zeros(length(x), 1);
    for i = 1: length(K)
        val(i) = Eur_Call_LVF_MC(S0, K(i), T, r, x, M, N) - V0_mkt(i);
    end
end

% Evaluate the objective function at x
F = Fx(x);

% Calculate the Jacobian only if needed
if nargout > 1
    col1 = ( Fx(x + (dx(1)*[1;0;0])) - F )/dx(1);
    col2 = ( Fx(x + (dx(2)*[0;1;0])) - F )/dx(2);
    col3 = ( Fx(x + (dx(3)*[0;0;1])) - F )/dx(3);
    J = [col1 col2 col3];
end

end