rel_err1_cov = [];
rel_err2_cov = [];
k1 = [];
k2 = [];
N = 10.*(2:8);
M = 10.*(2:8);
n = length(N);
m = length(M);
for i = 1:n
    for j = 1:m
        [rel_err1,rel_err2,k11,k22] = rel_err(N(i),M(j));
        rel_err1_cov = [rel_err1_cov rel_err1];
        rel_err2_cov = [rel_err2_cov rel_err2];
        k1 = [k1 k11];
        k2 = [k2 k22];
    end
end
rel_err1_cov = reshape(rel_err1_cov,m,n)';
rel_err2_cov = reshape(rel_err2_cov,m,n)';
figure(1)
surf(N,M,rel_err1_cov)
xlabel('N')
ylabel('M')
zlabel('relative error')
figure(2)
surf(N,M,rel_err2_cov)
xlabel('N')
ylabel('M')
zlabel('relative error')
figure(3)
plot(k1)
hold on
plot(k2)