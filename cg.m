function [u_cg,k] = cg(U0,f,A,tol,N,M)
U = [U0];
p0 = f-A*U;
r0 = p0;
r1 = r0'*r0;
k = 0;
while r1>tol
    Ap = A*p0;
    alpha = r1/(p0'*Ap);
    U = U+alpha*p0;
    r = r0-alpha*Ap;
    beta = r'*r/r1;
    p = r+beta*p0;
    r1 = r'*r;
    p0 = p;
    r0 = r;
    k = k+1;
end
u_cg = reshape(U,M-1,N-1);
U_cg = zeros(M+1,N+1);
U_cg(2:M,2:N) = u_cg;
u_cg = U_cg(:);