function [d, Q] = initialise_Hermite_transform_Golub_Welsch(N)
% Builds the orthogonal matrix Q and weight vector d such that
% the coeffs2vals transform is d .* (Q' * cfs) and
% the val2coeffs transform is Q * (vals ./ d).

% Calculate Q via eigendecomposition of Jacobi matrix:
beta = sqrt(.5*(1:N-1));
J = diag(beta, 1) + diag(beta, -1);
[Q, x] = eig(J,'vector');
[x, indx] = sort(x); % x = Gauss-Hermite quadrature nodes
Q = Q(:,indx);
% fix the sign of the columns of Q:
Q = Q .* (-1).^(mod(N,2) + (1:N)) .* sign(Q(N,:));

% d is proportional to abs(psi_{N-1}(x))
d = herm_func(N-1, x(floor(N/2)+1:end));
d = [d(end:-1:1); d((1+mod(N,2)):end)];
d = abs(d) * (sqrt(sum(d.^(-2) .* exp(-x.^2)))/pi^(1/4));
end

function val = herm_func(N, x)
% Evaluates the degree-N Hermite function up to an N-dependent constant
% Adapted from Chebfun's hermpts code due to Alex Townsend.

if N == 0
    val = exp(-x.^2/2);
elseif N == 1
    val = x.*exp(-x.^2/2);
elseif N <= 400 % evaluate using recurrence
    Hold = exp(-x.^2/2); H = sqrt(2)*x.*Hold;
    for k = 1:N-1
        val = x.*H.*sqrt(2/(k+1)) - Hold./sqrt(1+1/k);
        Hold = H; H = val;
    end
else % evaluate using Airy asymptotics DLMF (12.10.35)
    theta = acos(x./sqrt(2*N+1));
    musq = 2*N+1;
    cosT = cos(theta); sinT = sin(theta);
    sin2T = 2*cosT.*sinT;
    eta = .5*theta - .25*sin2T;
    chi = -(3*eta/2).^(2/3);
    phi = (-chi./sinT.^2).^(1/4);
    Airy0 = real(airy(musq.^(2/3)*chi));
    Airy1 = real(airy(1,musq.^(2/3)*chi));

    % Terms in DLMF (12.10.43):
    a0 = 1; b0 = 1;
    a1 = 15/144; b1 = -7/5*a1;
    a2 = 5*7*9*11/2/144^2; b2 = -13/11*a2;
    a3 = 7*9*11*13*15*17/6/144^3;

    % u polynomials in DLMF (12.10.9)
    u0 = 1; u1 = (cosT.^3-6*cosT)/24;
    u2 = (-9*cosT.^4 + 249*cosT.^2 + 145)/1152;
    u3 = (-4042*cosT.^9+18189*cosT.^7-28287*cosT.^5-151995*cosT.^3-259290*cosT)/414720;

    %first term
    A0 = 1;
    val = A0*Airy0;

    %second term
    B0 = -(a0*phi.^6.*u1+a1*u0)./chi.^2;
    val = val + B0.*Airy1./musq.^(4/3);

    %third term
    A1 = (b0*phi.^12.*u2 + b1*phi.^6.*u1 + b2*u0)./chi.^3;
    val = val + A1.*Airy0/musq.^2;

    %fourth term
    B1 = -(phi.^18.*u3 + a1*phi.^12.*u2 + a2*phi.^6.*u1 + a3*u0)./chi.^5;
    val = val + B1.*Airy1./musq.^(4/3+2);

    val = phi.*val;
end
end