%% A simple code to numerically solve Poisson's equation in 2D. 
close all; clear; clc
%Start with Dx = Dy

% Initiate with BCs, u = 0 on d\Omega, and \Omega (0,1)^2:
x_domain = [0,1]; y_domain = [0,1];
N = 100; % number of x gridpoints
M = 100; % number of y gridpoints
Dx = (x_domain(2)-x_domain(1))/N; Dy = (y_domain(2) - y_domain(1))/M;

% Make the actual Omega domain:
X_ = Dx:Dx:x_domain(2)-Dx; Y_ = Dy:Dy:y_domain(2)-Dy;
[X,Y] = ndgrid(X_,Y_);

% Define f, the inhomogenity of Possion's equation:
f = @(x_,y_) 13 * pi^2 * sin(2*pi*x_) .* sin(3*pi*y_);
f_grid = f(X,Y); f_vec = reshape(f_grid,[(N-1)*(M-1),1]);

% First need to code up the matrix for the scheme:
B = (1/Dy^2) * (diag(2*ones(1,M-1))+diag(-1*ones(1,M-2),1)+diag(-1*ones(1,M-2),-1)) + (2/Dx^2)*eye(M-1) ;
C = (1/Dx^2) * eye(M-1) ;


% Or straight from Wikipedia, Kronecker product.
nx = N-1; % number of grid points in the x-direction;
ny = M-1; % number of grid points in the y-direction;
ex = ones(nx,1);
Dxx = (1/Dx^2) * spdiags([ex -2*ex ex], [-1 0 1], nx, nx); %1D discrete Laplacian in the x-direction ;
ey = ones(ny,1);
Dyy = (1/Dy^2) * spdiags([ey, -2*ey ey], [-1 0 1], ny, ny); %1D discrete Laplacian in the y-direction ;
L = kron(Dyy, speye(nx)) + kron(speye(ny), Dxx) ;
L = -L; % Our Convention

% Backslash solver for U_vec and reshape onto grid.
U_vec = L\f_vec;
U_grid = reshape(U_vec,[M-1,N-1]);

% Whack the BCs back into the grid
U_grid = [zeros(1,(N-1)); U_grid; zeros(1,(N-1))];
U_grid = [zeros(M+1,1),U_grid,zeros(M+1,1)];

% Expand X,Y to include boundaries:
X_ = 0:Dx:x_domain(2); Y_ = 0:Dy:y_domain(2);
[X,Y] = ndgrid(X_,Y_);

% Plot the solution and observe,
figure(1)
surf(X,Y,U_grid); hold on

% Noting that the analytic solution is sin(2*pi*x)sin(3*pi*y) lets plot
% this too for comparison.
u_analytic = @(x,y) sin(2*pi*x) .* sin(3*pi*y);
U_Ana = u_analytic(X,Y);
figure(2)
surf(X,Y, U_Ana); hold on

figure(3)
error = U_Ana - U_grid;
surf(X,Y,error);
