function [u_ref_herm] = gross_pitaevskii_reference_sln(N,dt,T,beta,u0_fun)
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

filename=strcat('../data/u_ref_',dataset,'.mat');

if isfile(filename)
     load(filename,'u_ref_herm'); % File exists.
else
    %% Initialise Hermite transforms
    x = hermpts(N); % Only used for plotting purposes
    [d, Q] = initialise_Hermite_transform_Golub_Welsch(N);
    
    %% Choose initial condition
    % u0_fun = @(xx) exp(-0.5*(xx-1.0).^2) .* exp(1i*0.5*xx);
    
    u1_phys = u0_fun(x);
    
    % Convert to Hermite coefficients
    u1_herm = Q * (u1_phys./ d);
    
    %% Time evolution: Strang splitting
    
    for m = 1:M
        % ---- Half linear step in Hermite space ----
        u1_herm = exp(-i*dt*((0:N-1)'+1/2)) .* u1_herm;
    
        % ---- Full nonlinear step in physical space ----
    
        u1_phys = d .* (Q' * u1_herm);
        u1_phys = exp(-1i * dt * (beta * abs(u1_phys).^2)) .* u1_phys;
        u1_herm = Q * (u1_phys ./ d);
    
        % ---- Half linear kinetic step in Hermite space ----
        u1_herm = exp(-i*dt*((0:N-1)'+1/2)) .* u1_herm;
    end

    u_ref_herm=u1_herm;

    save(filename);
end
end