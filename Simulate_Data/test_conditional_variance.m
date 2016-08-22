

k = 2;
n = 1000;

rng(0);

x = randn( n,k);

eps = randn(n,1);

beta = rand(k,1);

y = x*beta + 3*eps;



cov(y)
conditional_variance_fn(y ,x)
