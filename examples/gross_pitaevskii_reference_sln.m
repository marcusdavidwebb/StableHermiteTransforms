function [psi_ref_herm] = gross_pitaevskii_reference_sln(N,dt,T,beta)
    % Returns reference solution of Gross-Pitaevskii equation based on
    % Strang splitting
    %
    % Ref solution loaded if already computed and stored in ../data/
    % Otherwise recompute with specified parameters
    % 
    % N...spatial discretisation param
    % dt...time step
    % T...end time
    % beta...nonlinearity strength

M = round(T/dt);
dt = T/M; % Correct for potential rounding discrepancy


%% Check if dataset already exists, if so load
dataset=strcat('N_',num2str(N),'_dt_',strrep(num2str(dt), '.', '-'),'_T_',strrep(num2str(T), '.', '-'),'_beta_',strrep(num2str(beta), '.', '-'));

filename=strcat('../data/psi_ref_',dataset,'.mat');

if isfile(filename)
     load(filename,'psi_ref_herm'); % File exists.
else
    %% Initialise Hermite transforms
    x = hermpts(N); % Only used for plotting purposes
    [d, Q] = initialise_Hermite_transform_Golub_Welsch(N);
    
    %% Choose initial condition
    psi0_fun = @(xx) exp(-0.5*(xx-1.0).^2) .* exp(1i*0.5*xx);
    
    psi1_phys = psi0_fun(x);
    
    % Convert to Hermite coefficients
    psi1_herm = Q * (psi1_phys./ d);
    
    %% Time evolution: Strang splitting
    
    for m = 1:M
        % ---- Half linear step in Hermite space ----
        psi1_herm = exp(-i*dt*((0:N-1)'+1/2)) .* psi1_herm;
    
        % ---- Full nonlinear step in physical space ----
    
        psi1_phys = d .* (Q' * psi1_herm);
        psi1_phys = exp(-1i * dt * (beta * abs(psi1_phys).^2)) .* psi1_phys;
        psi1_herm = Q * (psi1_phys ./ d);
    
        % ---- Half linear kinetic step in Hermite space ----
        psi1_herm = exp(-i*dt*((0:N-1)'+1/2)) .* psi1_herm;
    end

    psi_ref_herm=psi1_herm;

    save(filename);
end
end