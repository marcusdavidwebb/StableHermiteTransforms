function [T, Tinv] = initialise_Hermite_transform_unstable_v2(N)
%INITIALISE_HERMITE_TRANSFORM_DIRECT
% Build the Hermite transform matrices directly via the unstable
% three-term recurrence for Hermite functions:
%
%   psi_{j+1}(x) = sqrt(2/(j+1)) * x .* psi_j(x) - sqrt(j/(j+1)) * psi_{j-1}(x)
%
% Outputs:
%   T : NxN matrix, maps coefficients -> values
%               v = T * c
%
%   Tinv  : NxN matrix, maps values -> coefficients
%               c = Tinv * v
%
%   x         : Gauss-Hermite nodes
%   w         : Gauss-Hermite quadrature weights
%
% This is intended only as a baseline, since the recurrence is unstable
% for moderately large N.

    % Gauss-Hermite nodes and weights
    [x,w] = quad_gauss_hermite(N);

    % Preallocate matrix of Hermite function values:
    % Psi(k,m) = psi_{m-1}(x_k),   1 <= k,m <= N
    Psi = zeros(N, N);

    % Base cases:
    % psi_0(x) = pi^(-1/4) * exp(-x^2/2)
    Psi(:,1) = pi^(-1/4) * exp(-0.5 * x.^2);

    if N > 1
        % psi_1(x) = sqrt(2) * x * psi_0(x)
        Psi(:,2) = sqrt(2) * x .* Psi(:,1);
    end

    % Unstable recurrence for Hermite functions
    % Psi(:,m) corresponds to psi_{m-1}
    for m = 2:N-1
        Psi(:,m+1) = sqrt(2/m) * x .* Psi(:,m) ...
                   - sqrt((m-1)/m) * Psi(:,m-1);
    end

    % Backward transform:
    %   values = Hbackward * coeffs
    % so T(k,m) = psi_{m-1}(x_k)
    T = Psi;

    % Forward transform:
    % By Gauss-Hermite quadrature,
    %   c_m = sum_k w_k * exp(x_k^2) * psi_m(x_k) * f(x_k)
    % hence
    %   Hforward(m,k) = w(k) * exp(x(k)^2) * psi_{m-1}(x(k))
    %
    % Equivalently, row m is:
    %   (w .* exp(x.^2) .* Psi(:,m))'
    d = 1./(N*T(:,N).^2);
    Tinv = (Psi .* d).';
end