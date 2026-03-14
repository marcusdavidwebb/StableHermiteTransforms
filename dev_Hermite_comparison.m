clear all
rng(1) % Fix random seed for reproducability
%% Compare direct assembly of Hermite transform against stable one

addpath('modules/quadrature/')
addpath('modules/')

%% 
N=50;

[x,w]=quad_gauss_hermite(N); % quadrature points and nodes
[T,Tinv]=initialise_Hermite_transform_direct(N);
[d1,Q1]=initialise_Hermite_transform2(N);

norm(diag(d1)*Q1'-T)

norm(Q1*diag(1./d1)-Tinv)

norm(T*Tinv-eye(N))

% [d,Q]=initialise_Hermite_transform_old(x);
% a=hermpts(N)
% a=(a-a(N:-1:1))/2
% [d1,Q1]=initialise_Hermite_transform_old(a);
% 
% Q1
% norm(x-hermpts(N))


% 
% A=Q1'*Q1;
% 
% diag(A)
% A(end,end)