function [U1,U2]=uexact(N)
M=N;
a = 0; b = 1;
c = 0; d = 1;
% Define grid sizes
hx = (b-a)/(M); % length of sub-intervals in x-axis
hy = (d-c)/(N); % length of sub-intervals in y-axis
[X,Y] = meshgrid((a+hx):hx:(b-hx),(c+hy):hy:(d-hy));
uexact1 = @(X,Y)sin(2*pi*X).*sin(3*pi*Y);
uexact1=uexact1(X,Y);
U1=uexact1(:);

uexact2 = @(X,Y)(X-1).^5.*X.^2.*Y.*(Y-1);
uexact2=uexact2(X,Y);
U2=uexact2(:);

end