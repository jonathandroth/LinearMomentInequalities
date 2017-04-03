

k = 2;
n = 1000;

rng(0);

x = randn( n,k);
eps = randn(n,1);

beta = rand(k,1);

y = x*beta + 3*eps;

cov(y)
conditional_variance_fn(y ,x)
%%%%%%

%% Createa an equal matched pair on the x's

k = 2;
n = 1000;

rng(0);

x = randn( n,k) * [1,.5;.5,1]^(-1/2);
x = [x;x];
eps = randn(2* n,1);

beta = rand(k,1);

y = x*beta + 3*eps;

cov(y)
conditional_variance_fn(y ,x)

%%

N = 1000;
k = 5;


x = randn(N,k);
eps = randn(N,k);

y = x + eps;

unconditonal_var  = cov(y);

conditonal_var = conditional_variance_fn(y,x,0)
