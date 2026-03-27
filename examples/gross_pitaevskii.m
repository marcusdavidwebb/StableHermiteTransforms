%% Strang splitting for 1D Gross--Pitaevskii using Hermite transforms

clear all
clc

%% Paths
addpath('../modules/quadrature/')
addpath('../modules/')

%% Parameters

% Spatial / Hermite discretisation
N = 1024;

% Time discretisation
Tend = 5.0;
dt = 1e-3;%0.5*1e-4;
M = round(Tend/dt);
dt = Tend/M; % Correct for potential rounding discrepancy

% Gross--Pitaevskii parameters
beta = 1.0; % nonlinearity strength

dataset=strcat('N_',num2str(N),'_dt_',strrep(num2str(dt), '.', '-'),'_T_',strrep(num2str(Tend), '.', '-'),'_beta_',strrep(num2str(beta), '.', '-')); % Name of dataset for storage

%% Initialise Hermite transforms
x = hermpts(N); % Only used for plotting purposes
[T,Tinv] = initialise_Hermite_transform_unstable(N);
[d, Q] = initialise_Hermite_transform_Golub_Welsch(N);

%% Build differentiation / kinetic operator in Hermite space

% Example 1: shifted Gaussian with phase
alpha=0.125;
psi0_fun = @(xx) 1/sqrt(alpha)*exp(-alpha*(xx+25.0).^2) .* exp(1i*0.5*xx);

psi1_phys = psi0_fun(x);
psi2_phys = psi0_fun(x);

% Convert to Hermite coefficients
psi1_herm = Q * (psi1_phys./ d);
psi2_herm = Tinv * psi1_phys;


%% Time evolution: Strang splitting

for m = 1:M
    M-m
    % ---- Half linear step in Hermite space ----
    psi1_herm = exp(-i*dt*((0:N-1)'+1/2)) .* psi1_herm;

    % psi2_herm = exp(-i*dt*((0:N-1)'+1/2)) .* psi2_herm;

    % ---- Full nonlinear step in physical space ----

    psi1_phys = d .* (Q' * psi1_herm);
    psi1_phys = exp(-1i * dt * (beta * abs(psi1_phys).^2)) .* psi1_phys;
    psi1_herm = Q * (psi1_phys ./ d);

    % psi2_phys = T * psi2_herm;
    % psi2_phys = exp(-1i * dt * (beta * abs(psi2_phys).^2)) .* psi2_phys;
    % psi2_herm = Tinv * psi2_phys;

    % ---- Half linear kinetic step in Hermite space ----
    psi1_herm = exp(-i*dt*((0:N-1)'+1/2)) .* psi1_herm;

    % psi2_herm = exp(-i*dt*((0:N-1)'+1/2)) .* psi2_herm;
    
    if isnan(sum(psi1_herm))
        input('broke')
    end

    % Plot evolution of the solution field
    linestyle=["-","--","-.","-.","-.","-."];
    if mod(m-1,M/2-1)==0 %& m<3000
        psi1_phys = d .* (Q' * psi1_herm);
        % psi2_phys = d .* (Q' * psi2_herm);
        h(min(floor(m/(M/2-1))+1,3))=plot(x,abs(psi1_phys),linestyle(min(floor(m/(M/2-1))+1,3)),'Color','black','LineWidth',2,'MarkerFaceColor','white');
        hold on
        set(gca,'FontSize',16)
        xlabel('$x$','Interpreter','latex', 'FontSize', 22)
        ylabel('$|\psi(t)|$','Interpreter','latex', 'FontSize', 22)
        ylim([0,3])
        % xlim([-10,10])
        grid on
        
    end

    
end

% Add legend and save graph
legend(h,"t=0.0","t=2.5","t=5.0",'latex', 'FontSize', 16,'Location','northeast')
exportgraphics(gcf,strcat('../images/solution_field_GP_joint_',dataset,'.pdf'),'ContentType','vector')
