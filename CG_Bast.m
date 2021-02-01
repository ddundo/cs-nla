function[x,iter,X_values,residues] = CG_Bast(A, b, x, TOL)

    X_values = zeros(length(x),length(A));
    residues = zeros(1,length(A));
    
    res = b - A * x;
    p_vec = res;
    res_old = res' * res;
    iter = 0;
    for i = 1:length(b)
        iter = iter+1;
        Ap = A * p_vec;
        alpha = res_old / (p_vec' * Ap);
        x = x + alpha * p_vec;
        
        % To store an return the x_vectors computed on each iteration, for
        % analysis of convergence later
        X_values(:,i) = x;
        
        res = res - alpha * Ap;
        res_new = res' * res;
        residues(i) = res_new;
        if sqrt(res_new) < TOL
            X_values = X_values(:,1:iter);
              break
        end
        p_vec = res + (res_new / res_old) * p_vec;
        res_old = res_new;
    end
end