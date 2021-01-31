function [rel_err1,rel_err2,k1,k2] = rel_err(N,M)

%% Write the matrices
% N meshs in x direction, M meshs in y direction
% deltax = 1/N, deltay = 1/M
% £¨N = 100;£©
% £¨M = 100;£©
B = M^2*diag(2*ones(1,M-1))+M^2*diag((-1)*ones(1,M-2),1)+...
    M^2*diag((-1)*ones(1,M-2),-1)+2*N^2*eye(M-1);
C = N^2*eye(M-1);
tri_op = diag(ones(1,N-2),1)+diag(ones(1,N-2),-1);
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
u1_exact = [];
for i = 1:N+1
    for j = 1:M+1
        u1_exact = [u1_exact;fu1(x(i),y(j))];
    end
end
f2 = [];
for i = 2:N
    for j = 2:M
        f2 = [f2;ff2(x(i),y(j))];
    end
end
u2_exact = [];
for i = 1:N+1
    for j = 1:M+1
        u2_exact = [u2_exact;fu2(x(i),y(j))];
    end
end
%% Compute in backslash
u1_bs = A\f1;
u1_bs = reshape(u1_bs,M-1,N-1);
U1_bs = zeros(M+1,N+1);
U1_bs(2:M,2:N) = u1_bs;
u1_bs = U1_bs(:);

u2_bs = A\f2;
u2_bs = reshape(u2_bs,M-1,N-1);
U2_bs = zeros(M+1,N+1);
U2_bs(2:M,2:N) = u2_bs;
u2_bs = U2_bs(:);

%% Conjugate gradient method
tol = 10^(-4);
U0 = zeros((N-1)*(M-1),1);
[u1_cg,k1] = cg(U0,f1,A,tol,N,M);
[u2_cg,k2] = cg(U0,f2,A,tol,N,M);
%% error between conjugate gradient and backslash
rel_err1 = norm((u1_cg-u1_bs),Inf)/norm(u1_bs,Inf);
rel_err2 = norm((u2_cg-u2_bs),Inf)/norm(u2_bs,Inf);
end






