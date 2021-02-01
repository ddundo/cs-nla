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


% Define f, the inhomogenity of Possion's equation: THIS IS THE DEFINITION
% OF THE PROBLEM AT HAND!
f = @(x_,y_) -(x_-1).^3 .* (42*x_.^2 - 24.*x_ + 2) .* y_ .* (y_ - 1) - 2.*x_.^2.*(x_-1).^5; 
f_grid = f(X,Y); f_vec = reshape(f_grid,[(N-1)*(M-1),1]);
% We note the analytical solution for clarity here:
u_analytic = @(x,y) (x-1).^5 .* x.^2 .* y .* (y-1);
U_Ana = u_analytic(X,Y); U_Ana_vec = reshape(U_Ana, [(N-1)*(M-1),1]);


% First need to code up the matrix for the scheme:
B = (1/Dy^2) * (diag(2*ones(1,M-1))+diag(-1*ones(1,M-2),1)+diag(-1*ones(1,M-2),-1)) + (2/Dx^2)*eye(M-1) ;
C = (1/Dx^2) * eye(M-1) ;


% Or Kronecker product.
nx = N-1; % number of grid points in the x-direction;
ny = M-1; % number of grid points in the y-direction;
ex = ones(nx,1);
Dxx = (1/Dx^2) * spdiags([ex -2*ex ex], [-1 0 1], nx, nx); %1D discrete Laplacian in the x-direction ;
ey = ones(ny,1);
Dyy = (1/Dy^2) * spdiags([ey, -2*ey ey], [-1 0 1], ny, ny); %1D discrete Laplacian in the y-direction ;
L = kron(Dyy, speye(nx)) + kron(speye(ny), Dxx) ;
L = -L; % Our Convention


% CG solver for U_vec and reshape onto grid.
x0 = zeros((M-1)*(N-1),1); % Initial guess
[U_vec,no_iter,U_iters,RES_list] = CG_Bast(L,f_vec,x0,10^-5); % To Solve A*x = b, CG_Bast takes in inputs: CG_Bast(A,b,x0,TOL), where x_0 is the initial guess.
U_grid = reshape(U_vec,[M-1,N-1]);


%% Convergence analysis

% The theory indicates that |U - U^k| =< 2 |U-U^0| * (sqrt(k*A)) - 1) / (sqrt(k*A) + 1)
k_values = 1:no_iter; % each iteration
RES_list = RES_list(1:no_iter);

% The bound from theory:
kappa = condest(L);
c = log(2*sqrt((U_Ana_vec'*L*U_Ana_vec)));
m = log( (sqrt(kappa) - 1) / (sqrt(kappa) + 1) );
k_bound = m * k_values + c;

% The actual values:
A_Res = zeros(1,no_iter);
for i = 1:no_iter
    A_Res(i) = log(sqrt( (U_Ana_vec - U_iters(:,i))'*L*(U_Ana_vec - U_iters(:,i)) ));
end

log_RES = log(RES_list);

figure(1)
plot(k_values,k_bound,'linewidth',1.5); hold on;
plot(k_values,A_Res, 'linewidth',1.5); hold on;
legend("Theoretical Bound", "Actual A-norm residues");
title("Plot of convergence of CG method as bounded by the theoretical limit");
xlabel("K^{th} iterative step")
ylabel("log(||U - U^k||_A")
%semilogy(1:no_iter, RES_list(1:no_iter), 'linewidth', 1.5); hold on;
%% Now to plot
% Whack the BCs back into the grid
U_grid = [zeros(1,(N-1)); U_grid; zeros(1,(N-1))];
U_grid = [zeros(M+1,1),U_grid,zeros(M+1,1)];

% Expand X,Y to include boundaries:
X_ = 0:Dx:x_domain(2); Y_ = 0:Dy:y_domain(2);
[X,Y] = ndgrid(X_,Y_);
% Plot the solution and observe,
figure(2)
surf(X,Y,U_grid); hold on

% Noting that the analytic solution is sin(2*pi*x)sin(3*pi*y) lets plot
% this too for comparison - note that this is with expansion to include BCs
figure(3)
U_A = u_analytic(X,Y);
surf(X,Y,U_A)

% Error plotting
figure(4)
error = (U_A - U_grid)/norm(U_A);
surf(X,Y,error);
