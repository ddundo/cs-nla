function [u ,err, errvec]= GS2(A,f,u0,tol)
n = length(u0);
L = tril(A,-1);
R = triu(A,1);
D = spdiags(diag(A),0,n,n);
u = u0;
err = max(abs(A*u-f));
errvec =[err];
    while err > tol
    u = (D+L)\(-R*u+f);
    err = max(abs(A*u-f));
    errvec = [errvec err];
    end
end