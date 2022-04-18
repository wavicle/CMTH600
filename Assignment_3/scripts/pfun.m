function [F, J] = pfun(x)
% PFUN Returns the objective function value and the Jacobian
% Useage: [F, J] = pfun(x)
% x = input x vector of length 2
% F = input function value at x
% J = Jacobian of F at x

% First define F(x) as an inline function (row vector)
Fx = @(x) [
    3*(x(1)^2) + 2*x(1)*x(2) - 1; ...
    -3*x(2) + 5*x(1)*x(2) - 4; ...
    exp(x(1)) - sin(x(2)) + 1
];

% Define dx as a vector with two values (both = 0.01)
dx = [0.01, 0.01];

% Calculate the value for F(x) at given input x
F = Fx(x);

% Construct the Jacobian using two columns
col1 = ( Fx(x + (dx(1)*[1;0])) - Fx(x) )/dx(1);
col2 = ( Fx(x + (dx(2)*[0;1])) - Fx(x) )/dx(2);

J = [col1 col2];
end