function[x] = Jac_Bast(A, b, x0, TOL)
assert(length(x0) == length(b), "The length of b != x, Ax = b")
% Solve Ax = b, with initial guess x0 and 1000 iterations.

% JACOBI Iteration Matrix => Extract diagonal elements and leave remaining elements.
D = sparse(diag(diag(A))); RnL = sparse(A - D);

% We pre-compute the iteration matrix for efficiency
C_Jac = -D\(RnL);

% Lets see what the maximum lambda value is (useful for error stuff). Recall that the modulus needs to be below 1, 
maxLam = eigs(C_Jac, 1);
assert(abs(maxLam) < 1, "For the input matrix you have given, we will not obtain convergence, the spectral norm is larger than 1")

% Initialise the iteration.
x = x0;

count = 0;
res = 10; % Dummy initialisation
while res > TOL 
    x_past = x;
    x = C_Jac*x + D\b;
    res = norm((x - x_past), Inf);
    count = count + 1;
end

end