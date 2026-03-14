clear all
rng(1) % Fix random seed for reproducability
%% Compare direct assembly of Hermite transform against stable one

addpath('modules/quadrature/')
addpath('modules/')

%% 
N=11;

[x,w]=quad_gauss_hermite(N); % quadrature points and nodes
[T,Tinv]=initialise_Hermite_transform_direct(N);

% [d,Q]=initialise_Hermite_transform_old(x);
% a=hermpts(N)
% a=(a-a(N:-1:1))/2
% [d1,Q1]=initialise_Hermite_transform_old(a);
% 
% Q1
% norm(x-hermpts(N))


[d1,Q1]=initialise_Hermite_transform2(300);

A=Q1'*Q1;

diag(A)
% A(end,end)