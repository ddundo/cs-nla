function[x] = GS_Bast(A, b, x0, TOL)
assert(length(x0) == length(b), "The length of b != x, Ax = b")
% Solve Ax = b, with initial guess x0 and 1000 iterations.

% GAUSS-SEIDEL Iteration Matrix => Extract lower elements (D+L) for LHS, leave upper for RHS 
% Utilise sparse matrices for computational efficiency.
DnL = sparse(tril(A)); R = sparse(A - DnL); 

% We pre-compute the iteration matrix for efficiency
C_GS = -DnL\(R);

% Lets see what the maximum lambda value is (useful for error stuff). Recall that the modulus of this needs to be below 1
maxLam = eigs(C_GS, 1);
assert(abs(maxLam) < 1, "For the input matrix you have given, we will not obtain convergence, the spectral norm is larger than 1")
% Initialise the iteration.
x = x0;

count = 0;
res = 10; % Dummy initialisation
while res > TOL 
    x_past = x;
    x = C_GS*x + DnL\b;
    res = norm((x - x_past), Inf);
    count = count + 1;
end

end
