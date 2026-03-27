clear all;


%% Paths
addpath('../modules/quadrature/')
addpath('../modules/')

%psi0_fun = @(xx) exp(-0.5*(xx-1.0).^2) .* exp(1i*0.5*xx);
psi0_fun = @(xx) exp(-(xx-30.0).^2) .* 1./(3+cos(xx).^2);


N = 1024;

x = hermpts(N); % Only used for plotting purposes


[d, Q] = initialise_Hermite_transform_Golub_Welsch(N);


% psi0_fun = @(xx) exp(-(xx-35.0).^2) .* 1./(3+cos(xx).^2)%+exp(-(xx+20.0).^2) .* 1./(3+cos(xx).^2);


alpha=0.5;
psi0_fun = @(xx) 1/sqrt(alpha)*exp(-alpha*(xx+25.0).^2) .* exp(1i*0.5*xx);

psi1_phys = psi0_fun(x);

% Convert to Hermite coefficients
psi1_herm = Q * (psi1_phys./ d);

figure(1)
semilogy(abs(psi1_herm))

% figure(2)
% plot(x, real(psi1_phys))