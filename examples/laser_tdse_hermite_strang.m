clear all
clc

%% Set paths

addpath('../modules/quadrature/')
addpath('../modules/')

%% Set parameters
% Numerical method params
N = 128;
T = 200.0;
dt = 0.001;
M = round(T/dt);
store_every = 50;

% Physical params
Z = 1.0;
a = 1.0;
q = -1.0;

E0 = 0.02;
phi = 0.0;
omega = 0.25;
t_on = 0.0;
t_off = T;

%% Initialise Hermite transforms
[x, w] = quad_gauss_hermite(N);

[valweights, Q] = initialise_Hermite_transform(x);

A = (Q * (x .* Q')) - sqrt(2*(0:N-1)') .* diag(ones(N-1,1), -1); % Differentiation matrix
A(abs(A) < 1e-10) = 0;
A_s = sparse(A);

EXPM = expm(1i * (dt/2) * (A_s * A_s)); % Exponential \exp(i dt/2 d/dx^2)


%% Define initial condition

psi0_fun = @(xx) exp(-0.5*(xx-1.0).^2) .* exp(1i*0.5*xx);
psi_phys = psi0_fun(x);

norm0 = sqrt(sum(abs(psi_phys).^2 .* w));
psi_phys = psi_phys ./ norm0;
psi_herm = Q * (valweights .* psi_phys);


%% Define potentials
soft_core = @(xx) -Z ./ sqrt(xx.^2 + a.^2);

laser_E = @(t) (t > t_on && t < t_off) * (E0 * (sin(pi*(t - t_on)/(t_off - t_on))^2) * cos(omega*t + phi));
V_total = @(xx, t) soft_core(xx) - q * laser_E(t) * xx;

n_store = floor(M/store_every) + 1;
ts = zeros(n_store, 1);
norms = zeros(n_store, 1);
xexp = zeros(n_store, 1);
rho_store = zeros(n_store, N);

kstore = 1;
t = 0.0;
psi_phys = (valweights).^(-1) .* (Q' * psi_herm);

ts(kstore) = t;
norms(kstore) = sum(abs(psi_phys).^2 .* w);
xexp(kstore) = real(sum(conj(psi_phys) .* x .* psi_phys .* w));
rho_store(kstore, :) = abs(psi_phys).^2;


%% Evolution
for m = 1:M
    psi_herm = EXPM * psi_herm;

    t_mid = (m-0.5) * dt;
    psi_phys = (valweights).^(-1) .* (Q' * psi_herm);
    Vx = V_total(x, t_mid);
    psi_phys = exp(-1i * dt * Vx) .* psi_phys;
    psi_herm = Q * (valweights .* psi_phys);

    psi_herm = EXPM * psi_herm;
    norm(psi_herm)

    figure(1)
    plot(x,(valweights).^(-1).*(Q'*psi_herm))
    pause(0.1)





    % if mod(m, store_every) == 0
    %     t = m * dt;
    %     psi_phys = (valweights).^(-1) .* (Q' * psi_herm);
    %     kstore = kstore + 1;
    %     ts(kstore) = t;
    %     norms(kstore) = sum(abs(psi_phys).^2 .* w);
    %     xexp(kstore) = real(sum(conj(psi_phys) .* x .* psi_phys .* w));
    %     rho_store(kstore, :) = abs(psi_phys).^2;
    % end
end

figure(1); clf;
subplot(2,1,1);
plot(ts, norms, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('||\psi||^2');

subplot(2,1,2);
plot(ts, xexp, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('<x>');

figure(2); clf;
plot(x, rho_store(1,:), 'LineWidth', 1.5); hold on;
plot(x, rho_store(round(end/2),:), 'LineWidth', 1.5);
plot(x, rho_store(end,:), 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('|\psi(x,t)|^2');
legend({'t=0', sprintf('t=%.3g', ts(round(end/2))), sprintf('t=%.3g', ts(end))}, 'Location', 'best');

outdir = fullfile(repo_root, 'data');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
save(fullfile(outdir, 'laser_tdse_hermite_example.mat'), 'x', 'w', 'ts', 'norms', 'xexp', 'rho_store', 'N', 'T', 'dt', 'Z', 'a', 'q', 'E0', 'omega', 'phi', 't_on', 't_off');
