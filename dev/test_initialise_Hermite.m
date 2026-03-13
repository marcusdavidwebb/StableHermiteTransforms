
addpath('../modules/')
addpath('../modules/quadrature/')
n=10;
[x,w]=quad_gauss_hermite(n); % quadrature points and nodes


[d, Q] = initialise_Hermite_transform(x);



a=sqrt((1:(n-1))'/2)
A=spdiags(a,-1,n,n);
A=A+A';

[V,D]=eig(full(A));

V'+Q % Not zero but some columns match...

