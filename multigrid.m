%% Write the matrices
% N meshs in x direction, M meshs in y direction
% deltax = 1/N, deltay = 1/M
N = 32;
M = 32;
B = M^2*diag(2*ones(1,M-1))+M^2*diag((-1)*ones(1,M-2),1)+...
    M^2*diag((-1)*ones(1,M-2),-1)+2*N^2*eye(M-1);
C = N^2*eye(M-1);
tri_op = -diag(ones(1,N-2),1)-diag(ones(1,N-2),-1);
A = kron(eye(N-1),B)+kron(tri_op,C);
%% Example functions
x = linspace(0,1,N+1);
y = linspace(0,1,M+1);
ff1 = @(x,y) 13*pi^2*sin(2*pi*x).*sin(3*pi*y);
fu1 = @(x,y) sin(2*pi*x).*sin(3*pi*y); 
ff2 = @(x,y) -(x-1)^3*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5;
fu2 = @(x,y) (x-1).^5.*x.^2.*y.*(y-1);
f1 = [];
for i = 2:N
    for j = 2:M
        f1 = [f1;ff1(x(i),y(j))];
    end
end
f2 = [];
for i = 2:N
    for j = 2:M
        f2 = [f2;ff2(x(i),y(j))];
    end
end
%% Compute in backslash
u1_bs = A\f1;
u1_bs = reshape(u1_bs,M-1,N-1);
u1_bs = u1_bs(:);

u2_bs = A\f2;
u2_bs = reshape(u2_bs,M-1,N-1);
u2_bs = u2_bs(:);
%% 2-grid
L = tril(A,-1);
R = triu(A,1);
omega = 0.5;
tol = 10^(-6);
n = length(A);
D_inv = spdiags(1./diag(A),0,n,n);
    
u0 = zeros(length(A),1);
I1 = zeros((length(A)-1)/2,length(A));
for i = 1:(length(A)-1)/2
    I1(i,2*i-1:2*i+1) = [1 2 1];
end
I1 = 1/4*I1;
I2 = I1'*1/2;
A_bar = I1*A*I2;
err = [];

err1 = 1;
while err1 > tol
for i = 1:3
    u0 = ((1-omega)*eye(length(A))-omega*D_inv*(L+R))*u0+omega*D_inv*f1; 
end
r = f1-A*u0;
r_bar = I1*r;
e_bar = A_bar\r_bar;
e = I2*e_bar;
u0 = u0+e;
err1 = max(abs(A*u0-f1));
err = [err err1];
end

%% Multigrid
L = tril(A,-1);
R = triu(A,1);
omega = 0.5;
tol = 10^(-6);
n = length(A);
D_inv = spdiags(1./diag(A),0,n,n);

E = 1/8*eye(M-1)+1/16*diag(ones(M-2,1),1)+1/16*diag(ones(M-2,1),-1);
opt = zeros(N/2,M-1);
opt(1,1:2) = [2 1];
opt(N/2,end-1:end) = [1 2];
for i = 2:N/2-1
    opt(i,2*i-2:2*i) = [1 2 1];
end
I1 = kron(opt,E);
I2 = I1'*1/2;
A_bar = I1*A*I2;

% f1
u0 = zeros(length(A),1);
er1 = [];
err12 = [];
err1 = 1;
k12 = 0;
while err1 > tol
for i = 1:3
    u0 = ((1-omega)*eye(length(A))-omega*D_inv*(L+R))*u0+omega*D_inv*f1; 
end
r = f1-A*u0;
r_bar = I1*r;
e_bar = A_bar\r_bar;
e = I2*e_bar;
u0 = u0+e;
er1 = [er1 max(abs(u1_bs-u0))];
err1 = max(abs(A*u0-f1));
err12 = [err12 err1];
k12 = k12+1;
end

% f2
u0 = zeros(length(A),1);
er2 = [];
err22 = [];
err2 = 1;
k22 = 0;
while err2 > tol
for i = 1:3
    u0 = ((1-omega)*eye(length(A))-omega*D_inv*(L+R))*u0+omega*D_inv*f2; 
end
r = f2-A*u0;
r_bar = I1*r;
e_bar = A_bar\r_bar;
e = I2*e_bar;
u0 = u0+e;
er2 = [er2 max(abs(u2_bs-u0))];
err2 = max(abs(A*u0-f2));
err22 = [err22 err2];
k22 = k22+1;
end

figure(1)
semilogy(err12); hold on
semilogy(err22); 
title('convergence rate')
ylabel('|A*u_{mult}-f|')
legend('f1','f2')
xlabel('iterations')

figure(2)
semilogy(er1); hold on
semilogy(er2); 
title('convergence rate')
ylabel('|u_{exact}-u_{mult}|')
legend('f1','f2')
xlabel('iterations')


