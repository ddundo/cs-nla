%% A simple code to numerically solve Poisson's equation in 2D.
% This guy is a simple 2D solver for Poisson's equation in 2D, using a
% Simple Jacobi iteration method, for different M and N values.

close all; clear; clc

YSet = cell(10,10);
XSet = cell(10,10);
USet = cell(10,10);

h = figure;
filename = 'Poisson_Ex2_Variation.gif';
axis tight manual

%for N = 10:10:100 % number of x gridpoints
    %disp(N) % To indicate how far along we are.
    %for M = 10:10:100 % number of y gridpoints
        N = 50; M = 50;
        % Initiate with BCs, u = 0 on d\Omega, and \Omega (0,1)^2:
        x_domain = [0,1]; y_domain = [0,1];
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
        
        
        % Need to get the optimal value of w here.
        U_vec = GS_Bast(L,f_vec,x0,10^-8); % To Solve A*x = b
        
        
        U_grid = reshape(U_vec,[N-1,M-1]);
        
        % Expand U_grid to include BCs
        U_grid = [zeros(1,(M-1)); U_grid; zeros(1,(M-1))];
        U_grid = [zeros(N+1,1),U_grid,zeros(N+1,1)];
        
        % Expand X,Y to include boundaries:
        X_ = 0:Dx:x_domain(2); Y_ = 0:Dy:y_domain(2);
        [X,Y] = ndgrid(X_,Y_);
        
        U_A = u_analytic(X,Y);
        error = (U_A - U_grid)/norm(U_A);
        
        % Saving within cell arrays for animated GIF Later and error
        % analysis
        XSet{N/10,M/10} = X;
        YSet{N/10,M/10} = Y;
        USet{N/10,M/10} = error;
        
        % Draw as a GIF
        surf(X,Y,error)
        %drawnow
        % Capture frames as images
        %frame = getframe(h);
        %im = frame2im(frame);
        %[imind,cm] = rgb2ind(im,256);
        % Write GIF to file
        %if M/10 == 1 && N/10 == 1
        %    imwrite(imind,cm,filename,'gif','Loopcount',inf);
        %else
        %    imwrite(imind,cm,filename,'gif','WriteMode','append');
        %end
        
        figure(2)
        surf(X,Y,U_grid)
        
        figure(3)
        surf(X,Y,U_A)
       
    %end
%end

%msgbox("End of loop")