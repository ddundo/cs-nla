N=[32,64,128,256,512,1024];
size(N)

for ii=1:5
    tic
    [A, F]=Poisson(N(ii));
    b=F;
    U = zeros(size(F)); % first guess is zero vector
    r = b; % r = b - A*x starts equal to b
    p = r; % first search direction is r
    k=0;
    %E_0=sqrt((U-uexact(N(ii)))'*A*(U-uexact(N(ii))));
    while norm(r,2)>10^(-13)
        
        Ap=A*p; %construct and store
        alpha = r'*r / (p'*Ap);
        U = U + alpha * p; % step to next guess
        rnew = r - alpha * Ap; % update residual r
        beta = (rnew'*rnew)/(r'*r);
        r = rnew;
        p = r + beta * p; % compute new search direction
        k=k+1;
        ERROR(k)=norm(r,2);
        %ERROR(k)=sqrt((U-uexact(N(ii)))'*A*(U-uexact(N(ii))))/E_0;
    end
    k, toc
end

kappa=condest(A)
loglog(1:k,ERROR/ERROR(1)), hold on, loglog(1:k,((sqrt(kappa)-1)/(1+sqrt(kappa))).^(1:k))


