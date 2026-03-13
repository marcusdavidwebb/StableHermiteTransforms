function [d, Q] = initialise_Hermite_transform2(N)
    % Builds the orthogonal matrix Q and weight vector d such that
    % the coeffs2vals transform is d .* (Q' * cfs) and
    % the val2coeffs transform is Q * (vals ./ d).
    % Q is N x N and d is N x 1.
    
    % Calculate x and Q via eigendecomposition of Jacobi matrix:
    beta = sqrt(.5*(1:N-1));
    J = diag(beta, 1) + diag(beta, -1);
    [Q, D] = eig(J);
    [x, indx] = sort(diag(D));
     Q = Q(:,indx);
     Q = Q .* sign(Q(1,:)+eps);
     
     % Calculate d using asymptotics (valid for N > 200)
     
     if N >= 400
         theta0 = acos(x(floor(N/2)+1:end)./sqrt(2*N+1));
         [val, dval] = hermpoly_asy_airy(N, theta0);
         d = sqrt(2*N+1)*cos(theta0).*val + sqrt(2)*dval;
         d = [d(end:-1:1);d((1+mod(N,2)):end)];
     else
         [val, dval] = hermpoly_rec(N, x);
         d = x.*val + sqrt(2)*dval;
     end
     
     d = abs(d) * (sqrt(sum(d.^(-2) .* exp(-x.^2)))/pi^(1/4));
end

function [val, dval] = hermpoly_asy_airy(n, theta)
% HERMPOLY_ASY evaluation hermite poly using Airy asymptotic formula in
% theta-space.
% Taken from the internal Chebfun function due to Alex Townsend.

musq = 2*n+1;
cosT = cos(theta); sinT = sin(theta);
sin2T = 2*cosT.*sinT;
eta = .5*theta - .25*sin2T;
chi = -(3*eta/2).^(2/3);
phi = (-chi./sinT.^2).^(1/4);
const = 2*sqrt(pi)*musq^(1/6)*phi; 
Airy0 = real(airy(musq.^(2/3)*chi));
Airy1 = real(airy(1,musq.^(2/3)*chi));

% Terms in (12.10.43):
a0 = 1; b0 = 1;
a1 = 15/144; b1 = -7/5*a1;
a2 = 5*7*9*11/2/144^2; b2 = -13/11*a2;
a3 = 7*9*11*13*15*17/6/144^3;
b3 = -19/17*a3;

% u polynomials in (12.10.9)
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

val = const.*val;

%% Derivative

eta = .5*theta - .25*sin2T;
chi = -(3*eta/2).^(2/3);
phi = (-chi./sinT.^2).^(1/4);
const = sqrt(2*pi)*musq^(1/3)./phi;

% v polynomials in (12.10.10)
v0 = 1; v1 = (cosT.^3+6*cosT)/24;
v2 = (15*cosT.^4-327*cosT.^2-143)/1152;
v3 = (259290*cosT + 238425*cosT.^3 - 36387*cosT.^5 + 18189*cosT.^7 -...
    4042*cosT.^9)/414720;

%first term
C0 = -(b0*phi.^6.*v1 + b1.*v0)./chi;
dval = C0.*Airy0/musq.^(2/3);

% %second term
D0 =  a0*v0;
dval = dval + D0*Airy1;

% %third term
C1 = -(phi.^18.*v3 + b1*phi.^12.*v2 + b2*phi.^6.*v1 + b3*v0)./chi.^4;
dval = dval + C1.*Airy0/musq.^(2/3+2);

%fourth term
D1 = (a0*phi.^12.*v2 + a1*phi.^6.*v1 + a2*v0)./chi.^3;
dval = dval + D1.*Airy1/musq.^2;

dval = const.*dval;
end

function [val, dval] = hermpoly_rec(n, x0)
% HERMPOLY_rec evaluation of scaled Hermite poly using recurrence
x0 = x0*sqrt(2);
% evaluate:
Hold = exp(-x0.^2/4); H = x0.*exp(-x0.^2/4);
for k = 1:n-1
    Hnew = (x0.*H./sqrt(k+1) - Hold./sqrt(1+1/k));
    Hold = H; H = Hnew;
end
% evaluate derivative:
val = Hnew;
dval = (-x0.*Hnew + n^(1/2)*Hold);
end
