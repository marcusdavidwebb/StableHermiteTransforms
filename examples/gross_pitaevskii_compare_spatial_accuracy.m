%% Check performance of spatial discretisation with two types of Hermite transforms on the GP equation
clear all
clc

%% Paths
addpath('../modules/quadrature/')
addpath('../modules/')

%% Parameters
% Spatial / Hermite discretisation
N_ref=1024;

% Time discretisation
Tend = 5.0;
dt = 0.25*1e-3;
dt_ref = 0.25*1e-3;
M = round(Tend/dt);
dt = Tend/M; % Correct for potential rounding discrepancy

% Gross--Pitaevskii parameters
beta = 1.0;              % nonlinearity strength

% Initial condition: shifted Gaussian with phase
alpha=0.125;
u0_fun = @(xx) 1/sqrt(alpha)*exp(-alpha*(xx+25.0).^2) .* exp(1i*0.5*xx);
%% Compute/load reference solution
fprintf("Computing reference solution...\n")
u_ref_herm=gross_pitaevskii_reference_sln(N_ref,dt_ref,Tend,beta,u0_fun);

fprintf("Finished computing reference solution.\n")

%% Store solution accuracy
% N_vec=[4,6,8,11,16,23,32,45,64,91,128,181,256,362,512,724,1024,1448,2048,2896];
N_vec=20:20:1024;

%10:20:1024;

error_GW=zeros(size(N_vec));
error_direct=zeros(size(N_vec));

for jj=1:max(size(N_vec))
    N=N_vec(jj)

    %% Initialise Hermite transforms
    x = hermpts(N); % Only used for plotting purposes
    [T,Tinv] = initialise_Hermite_transform_unstable(N);
    [d, Q] = initialise_Hermite_transform_Golub_Welsch(N);
    
    %% Build differentiation / kinetic operator in Hermite space
    
    u1_phys = u0_fun(x);
    u2_phys = u0_fun(x);
    
    % Convert to Hermite coefficients
    u1_herm = Q * (u1_phys./ d);
    u2_herm = Tinv * u1_phys;
    
    
    %% Time evolution: Strang splitting
    
    for m = 1:M
        % ---- Half linear step in Hermite space ----
        u1_herm = exp(-i*dt*((0:N-1)'+1/2)) .* u1_herm;
    
        u2_herm = exp(-i*dt*((0:N-1)'+1/2)) .* u2_herm;
    
        % ---- Full nonlinear step in physical space ----
    
        u1_phys = d .* (Q' * u1_herm);
        u1_phys = exp(-1i * dt * (beta * abs(u1_phys).^2)) .* u1_phys;
        u1_herm = Q * (u1_phys ./ d);
    
        u2_phys = T * u2_herm;
        u2_phys = exp(-1i * dt * (beta * abs(u2_phys).^2)) .* u2_phys;
        u2_herm = Tinv * u2_phys;
    
        % ---- Half linear kinetic step in Hermite space ----
        u1_herm = exp(-i*dt*((0:N-1)'+1/2)) .* u1_herm;
    
        u2_herm = exp(-i*dt*((0:N-1)'+1/2)) .* u2_herm;
    end

    error_direct(jj) = sqrt(norm(u2_herm-u_ref_herm(1:N))^2+norm(u_ref_herm(N+1:end))^2);
    error_GW(jj) = sqrt(norm(u1_herm-u_ref_herm(1:N))^2+norm(u_ref_herm(N+1:end))^2);

end

%% Save data
dataset=strcat('N_',num2str(N),'_dt_',strrep(num2str(dt), '.', '-'),'_T_',strrep(num2str(Tend), '.', '-'),'_beta_',strrep(num2str(beta), '.', '-'));
filename=strcat('../data/u_ref_',dataset,'.mat');

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
ylabel('$L^2$-error at $t=5$','Interpreter','latex', 'FontSize', 22)
xlim([1,1000])
% ylim([1e-12,1])
legend(h,'Direct','Golub--Welsh','Interpreter','latex', 'FontSize', 16,'Location','northeast')
grid on
hold off

exportgraphics(gcf,strcat('../images/spatial_accuracy_GP_',dataset,'.pdf'),'ContentType','vector')
