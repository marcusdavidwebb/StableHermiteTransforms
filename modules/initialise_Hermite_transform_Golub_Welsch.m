function [d, Q] = initialise_Hermite_transform_Golub_Welsch(N)
%INITIALISE_HERMITE_TRANSFORM_GOLUB_WELSCH
% Builds the orthogonal matrix Q and weight vector d such that
% the coeffs2vals transform is d .* (Q' * cfs) and
% the val2coeffs transform is Q * (vals ./ d).
% See Webb--Maierhofer 2026.

J = diag(sqrt(.5:.5:(N-1)/2), 1) + diag(sqrt(.5:.5:(N-1)/2), -1);
[Q, x] = eig(J, 'vector');                % see Golub-Welsch 1969
[x, indx] = sort(x); Q = Q(:,indx);       % x = Gauss-Hermite nodes
Q = Q .* sign(Q(N,:)) .* (-1).^(N+1:2*N); % enforce signs of final row

d = sqrt(N) * abs(herm_func(N-1, x(floor(N/2)+1:end)));
d = [d(end:-1:1); d((1+mod(N,2)):end)];
end

function val = herm_func(N, x)
% Evaluates the degree-N Hermite function (for x in [0, sqrt(2N+1)])
if N <= 200 % use Clenshaw's algorithm
    val = exp(-x.^2/2) / pi^(1/4);
    val1 = zeros(size(x));
    for k = N:-1:1
        val2 = val1; val1 = val;
        val = x .* val1 * sqrt(2/k) - val2 / sqrt(1 + 1/k);
    end
else % use Airy asymptotics, DLMF (12.10.35), accurate for N > 200
    mu2 = 2*N+1; t = x/sqrt(mu2);      % 12.10.1
    theta = acos(t); t2 = t.^2;
    eta = (theta - t.*sqrt(1-t2))/2;   % 12.10.23
    zeta = -(3*eta/2).^(2/3);          % 12.10.39
    phi = (zeta./(t2-1)).^(1/4);       % 12.10.40

    % 12.10.43:
    a1 = 15/144; b1 = -7/5*a1;
    a2 = 5*7*9*11/2/144^2; b2 = -13/11*a2;
    a3 = 7*9*11*13*15*17/6/144^3;

    % 12.10.9:
    u1 = (t2-6).*t/24; u2 = ((-9*t2 + 249).*t2 + 145)/1152;
    u3 = ((((-4042*t2+18189).*t2-28287).*t2-151995).*t2-259290).*t/414720;

    % 12.10.42:
    phi6 = phi.^6; A0 = 1; B0 = -(phi6.*u1+a1)./zeta.^2;
    A1 = ((phi6.*u2 + b1*u1).*phi6 + b2)./zeta.^3;
    B1 = -(((phi6.*u3 + a1*u2).*phi6 + a2*u1).*phi6 + a3)./zeta.^5;

    % 12.10.35:
    Airy0 = airy(mu2^(2/3)*zeta);
    Airy1 = airy(1, mu2^(2/3)*zeta);
    val = Airy0.*(A0+A1/mu2^2) + (Airy1/mu2^(4/3)).*(B0+B1/mu2^2);

    % 12.10.14:
    g = (((-4027/4976640/mu2 + 1003/103680)/mu2 + 1/576)/mu2 - 1/24)/mu2 + 1;
    g = g*exp(-gammaln(N+1)/2+(mu2/4-1/12)*log(mu2)-(mu2-3)*log(2)/4-mu2/4);
    val = (pi^(1/4) * g) * phi .* val;
end
end