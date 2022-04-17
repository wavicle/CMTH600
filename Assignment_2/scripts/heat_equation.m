close all; clear all; clc;

M = 50; N = 50;
X = 0:1/M:1;
dtau = 1/N;

U = zeros(N,M + 1);
for i = 1: M
    U(1, i) = -(X(i)*X(i)) + X(i); % Initial condition
end
U(1 , M+1) = 0; % Boundary condition at x = 1

% Initialize some matrices for later use
I = eye(M+1);

vector1 = [1, dtau + X(1:M-1).^2, 1];
% vector2 and vector2 have an extra 0 each to support spdiags later below
vector2 = [0, 0, -dtau/2*ones(1,M-1)];
vector3 = [-dtau/2*ones(1,M-1), 0, 0];
% Use spdiags to create a sparse matrix, to save memory
M1 = spdiags(vector1', 0, M+1, M+1) ...
    + spdiags(vector2', 1, M+1, M+1) + spdiags(vector3', -1, M+1, M+1);

vector1 = [1, -dtau + X(1:M-1).^2, 1];
% vector2 and vector2 have an extra 0 each to support spdiags later below
vector2 = [0, 0, dtau/2*ones(1,M-1)];
vector3 = [dtau/2*ones(1,M-1), 0, 0];
% Use spdiags to create a sparse matrix, to save memory
M2 = spdiags(vector1', 0, M+1, M+1) ...
    + spdiags(vector2', 1, M+1, M+1) + spdiags(vector3', -1, M+1, M+1);

% Use LU decomposition for faster results
[L, U] = lu(M1);

for n = 1:N-1
    % Boundary conditions
    U(n+1, M) = U(n, M);
    % Solve using LU matrices calculated earlier for faster results
    B = M2*transpose(U(n, :));
    U(n+1, :) = transpose(U\(L\B));
end % end of N for-loop
