%% Write the matrices
% N meshs in x direction, M meshs in y direction
% deltax = 1/N, deltay = 1/M
N = 50;
M = 50;
B = M^2*diag(2*ones(1,M-1))+M^2*diag((-1)*ones(1,M-2),1)+...
    M^2*diag((-1)*ones(1,M-2),-1)+2*N^2*eye(M-1);
C_ssor = N^2*eye(M-1);
tri_op = diag(ones(1,N-2),1)+diag(ones(1,N-2),-1);
A = kron(eye(N-1),B)+kron(tri_op,C_ssor);
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
%% Splitting methods
D = diag(diag(A));
L = tril(A)- D;
R = triu(A)- D;
%% Jacobi
u_jc = zeros(length(A),1);
m = 500;
tol = 10^(-6);
k_jc = 1;
n = length(A);
while k_jc <= m
    err_jc = 0;
    u_j_0 = u_jc;
    for i = 1:n
        u_jc(i) = (-A(i,1:n)*u_j_0+f1(i)+A(i,i)*u_j_0(i))/A(i,i);
        err_jc = norm(u_jc-u_j_0,Inf);
    end  
    if err_jc <= tol
        break;
    else 
        k_jc = k_jc+1;
    end
end
%% Jacobi matrix version
u_j_1 = zeros(length(A),1);
m = 50;
tol = 10^(-6);
k_j_1 = 1;
E_j = D\eye(length(A));
C_j = -E_j*(R+L);
while k_j_1 < m
    u_j_0 = u_j_1;
    u_j_1 = C_j*u_j_0+E_j*f1; 
    err_jc_1 = norm(u_j_1-u_j_0,Inf);
    if err_jc_1 <= tol
        break;
    else 
        k_j_1 = k_j_1+1;
    end
end
%% Relaxed Jacobi matrix version
u_jr_1 = zeros(length(A),1);
m = 50;
tol = 10^(-6);
omega = 0.5;
E_jr = E_j;
C_jr = (1-omega)*eye(length(A))+omega*C_j;
while k_jr_1 < m
    u_jr_0 = u_jr_1;
    u_jr_1 = C_jr*u_jr_0+E_jr*f1; 
    err_jr_1 = norm(u_jr_1-u_jr_0,Inf);
    if err_jr_1 <= tol
        break;
    else 
        k_jr_1 = k_jr_1+1;
    end
end
%% Gauss-Seidel
% A = D+R+L
u_gs = zeros(length(A),1);
m = 500;
tol = 10^(-6);
k_gs = 1;
n = length(A);
while k_gs < m
    err_gs = 0;
    for i = 1:n
        diff = 0;
        diff = diff-A(i,1:n)*u_gs(1:n);
        diff = (diff+f1(i))/A(i,i);
        if abs(diff) > err_gs
            err_gs = abs(diff);
        end
        u_gs(i) = u_gs(i) +diff;
    end  
    if err_gs <= tol
        break;
    else 
        k_gs = k_gs+1;
    end
end
%% Gauss-Seidel matrix version
u_gs_1 = zeros(length(A),1);
m = 50;
tol = 10^(-6);
k_gs_1 = 1;
E_gs = (D+L)\eye(length(A));
C_gs = -E_gs*R;
while k_gs_1 < m
    u_gs_0 = u_gs_1;
    u_gs_1 = C_gs*u_gs_0+E_gs*f1; 
    err_gs_1 = norm(u_gs_1-u_gs_0,Inf);
    if err_gs_1 <= tol
        break;
    else 
        k_gs_1 = k_gs_1+1;
    end
end
%% SOR matrix version
u_sor_1 = zeros(length(A),1);
m = 100;
tol = 10^(-15);
Omega = linspace(1,2,50);
k = [];
for n = 1:length(Omega)
omega = Omega(n);
E_sor = (D+omega*L)\eye(length(A));
C_sor = E_sor*((1-omega)*D-omega*R);
E_sor = omega*E_sor;
k_sor_1 = 1;
while k_sor_1 < m
    u_sor_0 = u_sor_1;
    u_sor_1 = C_sor*u_sor_0+E_sor*f1; 
    err_sor_1 = norm(u_sor_1-u_sor_0,Inf);
    if err_sor_1 <= tol
        break;
    else 
        k_sor_1 = k_sor_1+1;
    end
end
k = [k k_sor_1];
end

%% SSOR matrix version
u_ssor_1 = zeros(length(A),1);
m = 300;
tol = 10^(-8);
omega = 1.5;
k_ssor_1 = 1;
C1 = (D+omega*R)\((1-omega)*D-omega*L);
C2 = (D+omega*L)\((1-omega)*D-omega*R);
C_ssor = C1*C2;
E_ssor = (D+omega*L)\eye(length(A))+(D+omega*R)\eye(length(A));
E_ssor = omega*E_ssor;
while k_sor_1 < m
    u_ssor_0 = u_ssor_1;
    u_ssor_1 = C_ssor*u_ssor_0+E_ssor*f1; 
    err_ssor_1 = norm(u_ssor_1-u_ssor_0,Inf);
    if err_ssor_1 <= tol
        break;
    else 
        k_ssor_1 = k_ssor_1+1;
    end
end