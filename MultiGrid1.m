global x_domain
x_domain = [0,1];
global y_domain
y_domain = [0,1];

% Example function 1 - the Laplacian
f1 = @(x,y) 13*pi^2 * sin(2*pi*x) .* sin(3*pi*y);
u1 = @(x,y) sin(2*pi*x) .* sin(3*pi*y);
% Example function 2 - the Laplacian
f2 = @(x,y) -(x-1).^3 .* (42*x.^2 - 24*x + 2) .* y.*(y-1) - 2*x.^2 .* (x-1).^5;
u2 = @(x,y) (x-1).^5 .* x.^2 .* y.*(y-1);

% Initialise with i = 5 - as in Demmel;
i = 5;
N = 2^i - 1; M = N;
Dx = (x_domain(2)-x_domain(1))/(N+1); Dy = (y_domain(2) - y_domain(1))/(M+1);
X_ = Dx:Dx:x_domain(2)-Dx; Y_ = Dy:Dy:y_domain(2)-Dy;

[X,Y] = ndgrid(X_,Y_);

% Example 1:
%f_grid = f1(X,Y); f_vec = reshape(f_grid,[(N)*(M),1]);
%x0 = zeros(length(f_vec),1);
%x = MGV(x0,f_vec,i);
%xGrid = reshape(x,[N,N]);
%figure(2)
%surf(X,Y,xGrid);

% Example 2
f_grid = f2(X,Y); f_vec = reshape(f_grid,[(N)*(M),1]);
x0 = zeros(length(f_vec),1);
x = MGV(x0,f_vec,i);
xGrid = reshape(x,[N,N]);
figure(3)
surf(X,Y,xGrid);


function[x_out] = MGV(x,b,i)
    disp(i)
    if i == 1
        % Work on this - should be able to solve exact problem.
        disp("Fill Me In")
        x_out = Solver(x,b,i);
    else
        disp(x)
        [x,T] = Solver(x,b,i);
        disp(x)
        r = T*x - b;
        
        % Plot error
        
        rPlot = reshape(r,[sqrt(length(r)),sqrt(length(r))]);
        rPlot = add_bnd(rPlot);
        [xPlot,yPlot] = meshgrid(linspace(0,1,sqrt(length(r)) + 2),linspace(0,1,sqrt(length(r)) + 2));
        figure(i+3)
        surf(xPlot,yPlot,rPlot); title("Residual at level " + num2str(i));
        
        rc = Restrict(r,i);
        d = Prolong(MGV(zeros(length(rc),1), rc, i-1),i-1); % 4.* not needed, though still out by a tad (1/20th or so)
        x = x - d;
        x_out = Solver(x,b,i);
    end

end


% Do a handful of relaxed Jacobi runs.
function[x_out, A] = Solver(x,b,i)
    
    A = make_A(2^i - 1, 2^i - 1);
    R = sparse(triu(A, 1)); 
    L = sparse(tril(A,-1));
    w = 3/4;
    Dinv = speye(size(A))./diag(A);
    LR = L + R;
    
    C = (1-w)*speye(size(A)) - w * Dinv *LR;
    
    f = w * Dinv *b;
    
    for j = 1:4
        x = C*x + f;
        
    end
    x_out = x;
    
end

% Restrict points to a coarse grid - naive, non-average.
function[RC] = Restrict(R,i)
    %assert(length(R) == (2^i - 1)^2);
    %disp(length(R))
    %disp((2^i - 1)^2)
    %disp(i)
    gridR = reshape(R, 2^i - 1, 2^i - 1);
    gridRC = gridR(2:2:end-1,2:2:end-1); % Naive implementation, pick out points.
    %size(gridRC)
    RC = reshape(gridRC, [(2^(i-1)-1)^2, 1]);
    %disp(RC)
end

% Prolong by interpolation - (not naive, interpolate)
function[E] = Prolong(EC, i)
    %assert(length(EC) == (2^i - 1)^2);
    gridEC = reshape(EC,[2^i-1, 2^i-1]);
    gridEC = add_bnd(gridEC);
    gridE = interp2(gridEC,1); %assert(length(gridE) == (2^(i+1)-1));
    gridE = gridE(2:end-1,2:end-1);
    E = reshape(gridE,[(2^(i+1)-1)^2, 1]);
end

function A = make_A(M, N)
    assert(M == N);
    global x_domain, global y_domain
    Dx = (x_domain(2)-x_domain(1))/(N+1); Dy = (y_domain(2) - y_domain(1))/(M+1);
    %disp(Dx)
    nx = N; % number of grid points in the x-direction;
    ny = M; % number of grid points in the y-direction;
    ex = ones(nx,1);
    Dxx = (1/Dx^2) * spdiags([ex -2*ex ex], [-1 0 1], nx, nx); %1D discrete Laplacian in the x-direction ;
    ey = ones(ny,1);
    Dyy = (1/Dy^2) * spdiags([ey, -2*ey ey], [-1 0 1], ny, ny); %1D discrete Laplacian in the y-direction ;
    A = kron(Dyy, speye(nx)) + kron(speye(ny), Dxx) ;
    A = -A; % Our Convention
    
end

function x = add_bnd(xx)
    x = zeros(size(xx) + 2);
    x(2:end-1, 2:end-1) = xx;
end
