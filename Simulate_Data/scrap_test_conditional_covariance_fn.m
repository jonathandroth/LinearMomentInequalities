
n= 1000;
k = 2;

%X = normrnd(0,1,n,k);
X = mvnrnd([0;0], [1,0.5;0.5,1], n);
Y = X + mvnrnd([0;0], [1,0.3;0.3,1], n);

cov(Y)
conditional_variance_fn(Y,X, 0)