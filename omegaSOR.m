    %% Best Omega for SOR
    w=0.1:0.05:1.99;

    iter1=0;iter2=0;
    for ii=1:length(w)
    [u, err, err_vecA] = SOR(w(ii), A, f, u0, tol);
    iter1(ii)=length(err_vecA);

    [u, err, err_vecB] = SOR(w(ii), A, g, u0, tol);
    iter2(ii)=length(err_vecB);
    end

    figure(3)
    subplot(1,2,2)
    semilogy(w,iter1,'r'),hold on,% semilogy(w,iter2,'b'), 
    %legend('SOR')
    xlabel('\omega'), ylabel('No of iteration')
    title('Varying of \omega on SOR iterations')

    
    
    %% Best Omega for SSOR
    w=0.1:0.05:1.99;

    iter5=0;iter6=0;
    for ii=1:length(w)
    [u, err, err_vecA] = SSOR(w(ii), A, f, u0, tol);
    iter5(ii)=length(err_vecA);

    [u, err, err_vecB] = SSOR(w(ii), A, g, u0, tol);
    iter6(ii)=length(err_vecB);
    end

    figure(3)
    subplot(1,2,2)
    semilogy(w,iter5,'k'),hold on, %semilogy(w,iter6,'b--')
    %legend('SSOR on Problem 1')
    xlabel('\omega'), ylabel('No of iteration')
    title('Varying of \omega on SSOR iterations')
