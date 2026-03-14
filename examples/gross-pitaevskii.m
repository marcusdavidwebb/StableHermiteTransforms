%% Strang splitting for 1D Gross--Pitaevskii using Hermite transforms

clear all
clc

rng(1)

%% Paths
addpath('../modules/quadrature/')
addpath('../modules/')

%% Parameters

% Spatial / Hermite discretisation
N = 50;

% Time discretisation
T = 20.0;
dt = 0.001;
M = round(T/dt);
store_every = 50;

% Gross--Pitaevskii parameters
beta = 1.0;              % nonlinearity strength
omega_trap = 0.5;        % harmonic trap frequency

% Optional time-dependent forcing:
use_forcing = false;
E0 = 0.05;
omega_laser = 0.5;
phi = 0.0;

%% Initialise Hermite transforms
x = hermpts(N); % Only used for plotting purposes
[T,Tinv] = initialise_Hermite_transform_unstable(N);
[d, Q] = initialise_Hermite_transform_Golub_Welsch(N);

%% Build differentiation / kinetic operator in Hermite space

% Example 1: shifted Gaussian with phase
psi0_fun = @(xx) exp(-0.5*(xx-1.0).^2) .* exp(1i*0.5*xx);

% Other possible examples:
% psi0_fun = @(xx) exp(-0.5*xx.^2);
% psi0_fun = @(xx) sech(xx);
% psi0_fun = @(xx) exp(-0.5*(xx+1).^2) + exp(-0.5*(xx-1).^2);

psi1_phys = psi0_fun(x);
psi2_phys = psi0_fun(x);

% Convert to Hermite coefficients
psi1_herm = Q * (psi1_phys./ d);
psi2_herm = Tinv * psi1_phys;


%% Time evolution: Strang splitting

for m = 1:M
    % ---- Half linear step in Hermite space ----
    psi1_herm = exp(-i*dt*((0:N-1)'+1/2)) .* psi1_herm;

    psi2_herm = exp(-i*dt*((0:N-1)'+1/2)) .* psi2_herm;

    % ---- Full nonlinear step in physical space ----

    psi1_phys = d .* (Q' * psi1_herm);
    psi1_phys = exp(-1i * dt * (beta * abs(psi1_phys).^2)) .* psi1_phys;
    psi1_herm = Q * (psi1_phys ./ d);

    psi2_phys = T * psi2_herm;
    psi2_phys = exp(-1i * dt * (beta * abs(psi2_phys).^2)) .* psi2_phys;
    psi2_herm = Tinv * psi2_phys;

    % ---- Half linear kinetic step in Hermite space ----
    psi1_herm = exp(-i*dt*((0:N-1)'+1/2)) .* psi1_herm;

    psi2_herm = exp(-i*dt*((0:N-1)'+1/2)) .* psi2_herm;

    
    %figure(1)
    if mod(m,100)==0
        psi1_phys = d .* (Q' * psi1_herm);
        psi2_phys = d .* (Q' * psi2_herm);
        plot(x,real(psi1_phys),'b-','LineWidth',2)
        hold on
        plot(x,real(psi2_phys),'r--','LineWidth',2)
        legend('Golub-Welsh','Direct')
        hold off
        ylim([-1,1])
        pause(0.05)
        m
    end
end