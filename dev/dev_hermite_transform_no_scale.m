clear all
rng(1) % Fix random seed for reproducability
%% Compare direct assembly of Hermite transform against stable one

addpath('modules/quadrature/')
addpath('modules/')

%% 
N=100;

%% Check if Gauss-Hermite weights are precomputed, if not use


filename=strcat('modules/quadrature/precomp/hermite_nodesweights',num2str(N),'.mat');

if isfile(filename)
     load(filename,'x','w'); % File exists.
else
    [x,w]=quad_gauss_hermite(N); % quadrature points and nodes
    save(filename,'x','w');    % File does not exist.
end

%%  

[d, Q] = initialise_Hermite_transform(x);

Hforward = Q'.*(1./d)';

%% Now compute directly
N = length(x);
Q = zeros(N);
hjm1 = ones(N,1) * pi^(-1/4);    % h_0(x) (first Hermite polynomial)
Q(:,1) = hjm1; 
if N > 1
    hj = sqrt(2) * x .* hjm1;    % h_1(x) (second Hermite polynomial)
    Q(:,2) = hj;
end
cum_log_scale = zeros(N);    % Initialize the cumulative log scaling matrix

% Loop to compute higher-order Hermite polynomials
for j = 3:N
    % Use recurrence relation for Hermite polynomials
    [hjm1, hj] = deal(hj, sqrt(2/(j-1)) * x .* hj - sqrt((j-2)/(j-1)) * hjm1);
    
    % Rescale values to avoid overflow
    %scale = arrayfun(@(v) (v < 100) * 1 + (v >= 100) * (1/abs(v)), abs(hj));
    %hj = hj .* scale; hjm1 = hjm1 .* scale;
    
    % Update cumulative log scaling
    %cum_log_scale(:,j) = cum_log_scale(:,j-1) + log(scale);

    % Assign the rescaled Hermite polynomial to the matrix
    Q(:,j) = hj;  
end

% Compute the weights
d = sqrt(N) * abs(Q(:,N)) .* exp(- x.^2 / 2);

% Normalize Q matrix (convert to unscaled Hermite functions)
for j = 1:N
    Q(:,j) = Q(:,j) ./ abs(Q(:,N)) .* exp(- 0.5 * log(N));
end

T=d.*Q;

norm(Hforward-Q^(-1).*(1./d)')
norm(Hforward-T^(-1))