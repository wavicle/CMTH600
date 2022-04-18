function [x, it, r] = Gauss_Newton(fun, x0, itmax, tol)
%
% Gauss-Newton method for solving nonlinear least squares problem without line search
%
% Input
% fun - function F(x), the objective function is 1/2*F'*F, which returns
% F(x) and corresponding Jacobian matrix J.
% x0 - initial value of x
% itmax - max number of iteration
% tol - stopping tolerance
%
% Output
% x - final result
% it - number of iterations
% r - objective function value 1/2*F'*F

r = []; % Values of the objective function for each iteration
it = 0; % Current iteration number
x = x0; % Start with initial vector = x0 and update it in each iteration

% We allow a maximum of itmax iterations
while it <= itmax
    [F, J] = fun(x);
    A = J'*J; g = J'*F;
    % Get descent direction by solving the equation A*h = -g
    h = A\-g;
    x = x + h; % Without line search, alpha = 1

    % Increase the iteration number
    it = it + 1;    
    % Note r for this iteration
    r(it) = 0.5*fun(x)'*fun(x);
    
    % Stopping criteria:
    % 1: Descent direction is close to 0
    % 2: Gradient is close to 0
    if norm(h) < tol || norm(g) < tol
        disp(['Terminating at iteration: ', num2str(it)]);
        break;
    end
end
