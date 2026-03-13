clear all
rng(1) % Fix random seed for reproducability
%% Compare direct assembly of Hermite transform against stable one

addpath('../modules/quadrature/')
addpath('../modules/')

%% 
N=100;

%% Check if Gauss-Hermite weights are precomputed, if not use


filename=strcat('../modules/quadrature/precomp/hermite_nodesweights',num2str(N),'.mat');

if isfile(filename)
     load(filename,'x','w'); % File exists.
else
    [x,w]=quad_gauss_hermite(N); % quadrature points and nodes
    save(filename,'x','w');    % File does not exist.
end

%%  

[d, Q] = initialise_Hermite_transform(x);
[d1, Q1] = initialise_Hermite_transform_unscaled(x);

norm(d-d1)
norm(Q-Q1)