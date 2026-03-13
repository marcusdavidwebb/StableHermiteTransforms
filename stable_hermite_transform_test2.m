clear all
rng(1) % Fix random seed for reproducability
%% Compare direct assembly of Hermite transform against stable one

addpath('modules/quadrature/')
addpath('modules/')


%% Initialise error measures
jj_max=90;
Nfactor=10;
time_new=zeros(jj_max,1);
time_direct=zeros(jj_max,1);
time_unscaled=zeros(jj_max,1);
error_matrixnorm=zeros(jj_max,1);
error_fn_approx_direct=zeros(jj_max,1);
error_fn_approx_new=zeros(jj_max,1);
error_fn_approx_unscaled=zeros(jj_max,1);
cond_unscaledforward=zeros(jj_max,1);
cond_Hforward=zeros(jj_max,1);
cond_Hbackward=zeros(jj_max,1);
cond_newforward=zeros(jj_max,1);
max_d_unscaled=zeros(jj_max,1);
min_d_unscaled=zeros(jj_max,1);
max_d=zeros(jj_max,1);
min_d=zeros(jj_max,1);
N_save=zeros(jj_max,1);


for jj=1:jj_max
    jj_max-jj
    %% Parameters N
    N=Nfactor*jj;
    N_save(jj)=N;

    %% Check if Gauss-Hermite weights are precomputed, if not use
    
    
    filename=strcat('modules/quadrature/precomp/hermite_nodesweights',num2str(N),'.mat');
    
    if isfile(filename)
         load(filename,'x','w'); % File exists.
    else
        [x,w]=quad_gauss_hermite(N); % quadrature points and nodes
        save(filename,'x','w');    % File does not exist.
    end
    
    %% Unscaled Hermite transform
    tic
    [d1, Q1] = initialise_Hermite_transform_unscaled(x);
    time_unscaled(jj)=toc;
    
    %% New Hermite transform
    tic
    [d, Q] = initialise_Hermite_transform(x);
    time_new(jj)=toc;

    %% Some tests
    
    % Error in matrix norm

    error_matrixnorm(jj)=norm(Q'.*(1./d)' - Q1'.*(1./d1)');
    

    %norm(Q1.*valweights' - Hforward)
    
    % Error in function approximation and return
    f=@(x) 1./(2+cos(x)).*exp(-x.^2/2);
    %f=@(x) sum(exp(i*x*rand(1,3)),2).*exp(-x.^2/2);
    f_val=f(x);
    
    f_coeff1=Q1'*((1./d1).*f_val);
    f_coeff2=Q'*((1./d).*f_val);
    
    f_val1=d1.*Q1*(f_coeff2);
    f_val2=d.*Q*(f_coeff2);
    
    error_fn_approx_unscaled(jj)=norm(f_val-f_val1);
    error_fn_approx_new(jj)=norm(f_val-f_val2);

    cond_newforward(jj)=cond(Q'.*(1./d)');
    cond_unscaledforward(jj)=cond(Q1'.*(1./d1)');
    % cond_Hbackward(jj)=cond(Hbackward);
    max_d(jj)=min(abs(d));
    min_d(jj)=max(abs(d));
    max_d_unscaled(jj)=min(abs(d1));
    min_d_unscaled(jj)=max(abs(d1));
    
    % To be continued
    save(strcat("data/numerical_eval_N_",num2str(Nfactor),"-",num2str(jj_max*Nfactor),".mat"))

end