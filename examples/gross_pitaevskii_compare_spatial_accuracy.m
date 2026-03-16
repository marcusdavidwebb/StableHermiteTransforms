%% Check performance of spatial discretisation with two types of Hermite transforms on the GP equation
clear all
clc

%% Paths
addpath('../modules/quadrature/')
addpath('../modules/')

%% Parameters
% Spatial / Hermite discretisation
N_ref=256;

% Time discretisation
Tend = 1.0;
dt = 1e-3;
dt_ref = 1e-3;
M = round(Tend/dt);
dt = Tend/M; % Correct for potential rounding discrepancy

% Gross--Pitaevskii parameters
beta = 1.0;              % nonlinearity strength

% Initial condition: shifted Gaussian with phase
psi0_fun = @(xx) exp(-0.5*(xx-1.0).^2) .* exp(1i*0.5*xx);

%% Compute/load reference solution
fprintf("Computing reference solution...\n")
psi_ref_herm=gross_pitaevskii_reference_sln(N_ref,dt_ref,Tend,beta);

fprintf("Finished computing reference solution.\n")

%% Store solution accuracy
% N_vec=[4,6,8,11,16,23,32,45,64,91,128,181,256,362,512,724,1024,1448,2048,2896];
N_vec=4:1:100;

error_GW=zeros(size(N_vec));
error_direct=zeros(size(N_vec));

for jj=1:97
    N=N_vec(jj)

    %% Initialise Hermite transforms
    x = hermpts(N); % Only used for plotting purposes
    [T,Tinv] = initialise_Hermite_transform_unstable(N);
    [d, Q] = initialise_Hermite_transform_Golub_Welsch(N);
    
    %% Build differentiation / kinetic operator in Hermite space
    
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
    end

    error_direct(jj) = sqrt(norm(psi2_herm-psi_ref_herm(1:N))^2+norm(psi_ref_herm(N+1:end))^2);
    error_GW(jj) = sqrt(norm(psi1_herm-psi_ref_herm(1:N))^2+norm(psi_ref_herm(N+1:end))^2);

end

%% Save data
dataset=strcat('N_',num2str(N),'_dt_',strrep(num2str(dt), '.', '-'),'_T_',strrep(num2str(Tend), '.', '-'),'_beta_',strrep(num2str(beta), '.', '-'));
filename=strcat('../data/psi_ref_',dataset,'.mat');

save(filename);

%% Plotting
clear h
figure(1)
h(2)=semilogy(N_vec,error_GW,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(1)=semilogy(N_vec,error_direct,'--','Color','#c74cb9','LineWidth',3,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
legend('Direct','Golub-Welsh')


set(gca,'FontSize',16)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)
ylabel('$L^2$-error at $t=1$','Interpreter','latex', 'FontSize', 22)
xlim([1,100])
ylim([1e-12,1])
legend(h,'Direct','Golub--Welsh','Interpreter','latex', 'FontSize', 16,'Location','northeast')
grid on
hold off

exportgraphics(gcf,strcat('../images/spatial_accuracy_GP_',dataset,'.pdf'),'ContentType','vector')
