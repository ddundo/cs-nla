function [A, F]=Poisson(N)
M = N; % number of points
% Define domain
a = 0; b = 1;
c = 0; d = 1;
% Define grid sizes
hx = (b-a)/(M); % length of sub-intervals in x-axis
hy = (d-c)/(N);

%% Poisson matrix
r2 = 2*ones(N-1,1);
r = -ones(N-2,1);
B = diag(r2,0) + diag(r,1) + diag(r,-1);
% Sparse matrix B
B = sparse(B);full(B);
% Build sparse identity matrix
I = speye(N-1);full(I);
% Build tridiagonal block matrix A
A = hx^(-2)*(kron(B,I) + kron(I,B));

%% f without BC
% Generate 2D arrays of grids
[X,Y] = meshgrid((a+hx):hx:(b-hx),(c+hy):hy:(d-hy));
f = @(X,Y)13*pi^2*sin(2*pi*X).*sin(3*pi*Y);
F=f(X(:), Y(:));
end