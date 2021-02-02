function Ue=uexact(N)
M=N;
a = 0; b = 1;
c = 0; d = 1;
% Define grid sizes
hx = (b-a)/(M); % length of sub-intervals in x-axis
hy = (d-c)/(N); % length of sub-intervals in y-axis
[X,Y] = meshgrid((a+hx):hx:(b-hx),(c+hy):hy:(d-hy));
uexact = @(X,Y)sin(2*pi*X).*sin(3*pi*Y);
uexact=uexact(X,Y);
Ue=uexact(:);
end