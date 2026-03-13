%% Code that computes scattering operator Su_- for a given u_-

clear all
rng(1) % Fix random seed for reproducability
addpath('modules/quadrature/')
addpath('modules/gauss_legendre_quad/')
addpath('modules/')
addpath('modules/Hermite_transforms/')

%% Input function psi

psi0=@(x) exp(i*x).*exp(-(x-1.0).^2/2)+exp(-(x+2.0).^2/4); % Stable solution??

%% Parameters N
N_ref=1024;
N=512;
deltaa=1.0;
T=0.1;
mu=1.0;

%% Potential term
%   We solve:
%   psi_t = i psi_xx - i V(psi,x) psi,
potential_term="cubicNLSE"; %, "quartic_potential", "quadratic_potential", "double_well_potential"

if potential_term=="cubicNLSE"
    V=@(psi,x) -mu*abs(psi).^2;
elseif potential_term=="quartic_potential"
    V=@(psi,x) 1/8*abs(x).^4;
elseif potential_term=="quadratic_potential"
    V=@(psi,x) abs(x).^2;
elseif potential_term=="double_well_potential"
    a=1/4;
    b=1;
    V=@(psi,x) a.^2.*(abs(x).^2-b.^2).^2;
end

%% Setting up the Hermite transforms


filename_ref=strcat('modules/quadrature/precomp/hermite_nodesweights',num2str(N_ref),'.mat');

if isfile(filename_ref)
     load(filename_ref,'x','w'); % File exists.
else
    [x,w]=quad_gauss_hermite(N_ref); % quadrature points and nodes
    save(filename_ref,'x','w');    % File does not exist.
end

x_ref=x;
w_ref=w;

filename=strcat('modules/quadrature/precomp/hermite_nodesweights',num2str(N),'.mat');

if isfile(filename)
     load(filename,'x','w'); % File exists.
else
    [x,w]=quad_gauss_hermite(N); % quadrature points and nodes
    size(x)
    save(filename,'x','w');    % File does not exist.
end

% Stable transform:

[valweights_ref, Q_ref] = initialise_Hermite_transform(N_ref,x_ref);
[valweights, Q] = initialise_Hermite_transform(N,x);

% From fn vals to coeffs: Q*(valweights.*uminus_phys)
% From coeffs to fn vals: (valweights).^(-1).*(Q'*u_coeff)

%% Set up sparse differentiation matrix

A=(Q*(x.*Q'))-sqrt(2*(0:N-1)').*diag(ones(N-1,1),-1);

A(abs(A)<1e-10)=0;

A_s=sparse(A);

% For reference solution:


A_ref=(Q_ref*(x_ref.*Q_ref'))-sqrt(2*(0:N_ref-1)').*diag(ones(N_ref-1,1),-1);

A_ref(abs(A_ref)<1e-10)=0;

A_s_ref=sparse(A_ref);

% From f to df/dx is multiplication with A_s

%% Set up function

psi0_phys=psi0(x);
psi0_ref_phys=psi0(x_ref);
psi0_herm=Q*(valweights.*psi0_phys); 
psi0_ref_herm=Q_ref*(valweights_ref.*psi0_ref_phys); 


%% Hermite interpolation

H=@(n,x) exp(-1/2*x.^2).*hermiteH(n,x)./sqrt(2.^n.*sqrt(pi).*factorial(n));

% herm_interpolated = @(x) transpose(psi0_herm)*H(1:N,x);
% 
% function output=herm_interpolated(N,x)
%     for j=0:N-1
%         output(j) = herm_interpolated(j,x(j));
%     end
% end

sqrt(integral(@(x) abs(herm_interpolated(N,psi0_herm,x)).^2,-20,20))

%% Reference solution
M_ref=T*1e4;%2*1e6;
tau_ref=T/M_ref;


ref_filename=strcat('data/ref_solution_',potential_term,'_N_',num2str(N_ref),'_T_',replace(num2str(T),'.','-'),'_tau_',replace(num2str(tau_ref),'.','-'),'.mat');

if isfile(ref_filename)
    load(ref_filename,'psi_ref_herm');
else
    psi_ref_herm=psi0_ref_herm;
    
    EXPM_ref=expm(i*tau_ref/2*A_s_ref*A_s_ref);
    for m=1:M_ref
        M_ref-m
        psi_ref_herm=EXPM_ref*psi_ref_herm;
        
        
        psi_ref_phys=(valweights_ref).^(-1).*(Q_ref'*psi_ref_herm);
        
        aux_ref=V(psi_ref_phys,x_ref);
        
        psi_ref_phys=exp(-i*tau_ref*aux_ref).*psi_ref_phys;
        
        psi_ref_herm=Q_ref*(valweights_ref.*psi_ref_phys);


        psi_ref_herm=EXPM_ref*psi_ref_herm;
    end
    
    save(ref_filename,'psi_ref_herm');
end


%% Estimate convergence rates

global_error_lie=zeros(1,10);
global_error_strang=zeros(1,10);
time_lie=zeros(1,10);
time_strang=zeros(1,10);
tau_jj=zeros(1,10);


for jj=1:10
    M=2^jj;
    jj

    %tau=0.0001;
    %M=floor(T/tau);
    tau=T/M;
    %M=1;
    % M=100
    %% And next...integrate the equation
    
    psi_herm=psi0_herm;
    u_herm=psi0_herm;
    
    % Linear part
    EXPM=expm(i*tau/2*A_s*A_s);

    tic
    for m=1:M
        M-m
        % Lie splitting
        u_herm=EXPM*EXPM*u_herm;
        u_phys=(valweights).^(-1).*(Q'*u_herm);
        aux_u=V(u_phys,x);
        u_phys=exp(-i*tau*aux_u).*u_phys;
        u_herm=Q*(valweights.*u_phys);
    end
    time_lie(jj)=toc;

    tic
    for m=1:M
        M-m
        % Strang splitting
        psi_herm=EXPM*psi_herm;
        
        
        psi_phys=(valweights).^(-1).*(Q'*psi_herm);
        
        aux=V(psi_phys,x);%mu*psi_phys.*conj(psi_phys);
        
        psi_phys=exp(-i*tau*aux).*psi_phys;
        
        psi_herm=Q*(valweights.*psi_phys);


        psi_herm=EXPM*psi_herm;
    end
    
    time_strang(jj)=toc;
    global_error_strang(jj)=sqrt(norm(psi_ref_herm(1:N)-psi_herm).^2+norm(psi_ref_herm(N+1:end)).^2);
    global_error_lie(jj)=sqrt(norm(psi_ref_herm(1:N)-u_herm).^2+norm(psi_ref_herm(N+1:end)).^2);
    %% Now try direct approach
    


    %% Check consistency
    
    %global_error(jj)=norm(Q*(valweights.*(psi1_phys_new-psi11_phys)));
    tau_jj(jj)=tau;
end

%% Save error data
data_name=strcat('data/Hermite_error_',potential_term,'_Nref_',num2str(N_ref),'_N_',num2str(N),'_T_',replace(num2str(T),'.','-'),'_tauref_',replace(num2str(tau_ref),'.','-'),'.mat');
save(data_name,'global_error_strang','global_error_lie','time_strang','time_lie','tau_jj');
%% Plot this
figure(2)
loglog(tau_jj,global_error_lie, 'LineWidth',2)
hold on
loglog(tau_jj,global_error_strang, 'LineWidth',2)
grid on
loglog(tau_jj,tau_jj,'r--', 'LineWidth',2)
loglog(tau_jj,20*tau_jj.^2,'r--', 'LineWidth',2)
xlabel('$\tau$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(['$\|\psi_{\mathrm{direct}}-\psi_{R}\|_{L^2}$'], 'Interpreter', 'latex', 'FontSize', 16);
legend({'Lie','Strang' '$\mathcal{O}(\tau,\tau^2)$'}, ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
grid on;
ylim([1e-9,10])
hold off;

% Save the figure if desired
%exportgraphics(gcf, strcat('images/global_error_',potential_term,'_T_',replace(num2str(T),'.','-'),'.pdf'), 'ContentType', 'vector');

%% CPU time vs error

figure(3)
loglog(time_lie,global_error_lie, 'LineWidth',2)
hold on
loglog(time_strang,global_error_strang, 'LineWidth',2)
grid on
xlabel('CPU-time', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(['$\|\psi_{\mathrm{direct}}-\psi_{R}\|_{L^2}$'], 'Interpreter', 'latex', 'FontSize', 16);
legend({'CN','new'}, ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
grid on;
ylim([1e-9,10])
hold off;

% Save the figure if desired
%exportgraphics(gcf, strcat('images/global_error_vs_CPU_time_',potential_term,'_T_',replace(num2str(T),'.','-'),'.pdf'), 'ContentType', 'vector');

%% Final solution value plot

% Reconstruct final physical solutions on their native Hermite nodes:
u_phys_final      = (valweights).^(-1)     .* (Q'     * u_herm);        % Lie final on x
psi_phys_final    = (valweights).^(-1)     .* (Q'     * psi_herm);      % Strang final on x
psi_ref_phys_full = (valweights_ref).^(-1) .* (Q_ref' * psi_ref_herm);  % Ref final on x_ref

% Interpolate reference onto coarse Hermite nodes x (so curves are comparable)
[x_ref_sorted, idxR] = sort(x_ref);
psi_ref_sorted       = psi_ref_phys_full(idxR);
psi_ref_on_x         = interp1(x_ref_sorted, psi_ref_sorted, x, 'spline');

% Guard against NaNs (should be rare for Hermite nodes, but safe)
psi_ref_on_x(isnan(psi_ref_on_x)) = 0;

figure(4); clf;

% ---- Plot magnitude |psi| ----
subplot(2,1,1);
plot(x, abs(psi0_phys), 'LineWidth', 1.5); hold on;          % initial on x
plot(x_ref, abs(psi_ref_phys_full), 'LineWidth', 2.0);       % ref on x_ref
plot(x, abs(psi_ref_on_x), ':', 'LineWidth', 2.0);           % ref interpolated to x
plot(x, abs(u_phys_final), '--', 'LineWidth', 1.8);          % Lie final on x
plot(x, abs(psi_phys_final), '--', 'LineWidth', 1.8);        % Strang final on x
grid on;
xlabel('$x$', 'Interpreter','latex', 'FontSize', 14);
ylabel('$|\psi(x)|$', 'Interpreter','latex', 'FontSize', 14);
title(sprintf('Final solutions at T = %.3g (N=%d, N_{ref}=%d, finest \\tau=%.3g)', ...
      T, N, N_ref, tau_jj(end)), 'Interpreter','latex', 'FontSize', 14);

legend({ '$|\psi_0|$ on $x$', ...
         '$|\psi_{\mathrm{ref}}(T)|$ on $x_{\mathrm{ref}}$', ...
         '$|\psi_{\mathrm{ref}}(T)|$ interpolated to $x$', ...
         '$|\psi_{\mathrm{Lie}}(T)|$', ...
         '$|\psi_{\mathrm{Strang}}(T)|$' }, ...
       'Interpreter','latex', 'FontSize', 11, 'Location', 'best');

% ---- Plot pointwise error vs interpolated reference on x ----
subplot(2,1,2);
plot(x, abs(u_phys_final   - psi_ref_on_x), 'LineWidth', 1.8); hold on;
plot(x, abs(psi_phys_final - psi_ref_on_x), 'LineWidth', 1.8);
grid on;
xlabel('$x$', 'Interpreter','latex', 'FontSize', 14);
ylabel('$|\psi(T)-\psi_{\mathrm{ref}}(T)|$', 'Interpreter','latex', 'FontSize', 14);
legend({'Lie error on $x$', 'Strang error on $x$'}, ...
       'Interpreter','latex', 'FontSize', 12, 'Location', 'best');

% Optional export
% exportgraphics(gcf, strcat('images/final_solution_',potential_term,'_T_',replace(num2str(T),'.','-'),'.pdf'), ...
%               'ContentType', 'vector', 'Bounds', 'loose');
