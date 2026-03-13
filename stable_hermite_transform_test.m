clear all
rng(1) % Fix random seed for reproducability
%% Compare direct assembly of Hermite transform against stable one

addpath('modules/quadrature/')
addpath('modules/')


%% Initialise error measures
jj_max=15;
Nfactor=10;
time_new=zeros(jj_max,1);
time_direct=zeros(jj_max,1);
error_matrixnorm=zeros(jj_max,1);
error_fn_approx_direct=zeros(jj_max,1);
error_fn_approx_new=zeros(jj_max,1);
cond_Hforward=zeros(jj_max,1);
cond_Hbackward=zeros(jj_max,1);
cond_newforward=zeros(jj_max,1);
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
    
    
    %% Setting up the Hermite transforms in classical form
    
    H=@(n,x) exp(-1/2*x.^2).*hermiteH(n,x)./sqrt(2.^n.*sqrt(pi).*factorial(n));
    H_poly=@(n,x) hermiteH(n,x)./sqrt(2.^n.*sqrt(pi).*factorial(n));
    H_half=@(n,x) exp(-1/4*x.^2).*hermiteH(n,x)./sqrt(2.^n.*sqrt(pi).*factorial(n));
    
    % Forward transform
    
    Hforward=zeros(N);
    tic
    for m=1:N
        for k=1:N
            Hforward(m,k)=H_poly(m-1,x(k))*w(k)*exp(x(k)^2/2);
        end
    end
    time_direct(jj)=toc;
    % Inverse transform
    
    Hbackward=zeros(N);
    for k=1:N
        for m=1:N
            Hbackward(k,m)=H(m-1,x(k));
        end
    end
    
    
    %% New Hermite transforms
    tic
    [d, Q] = initialise_Hermite_transform(x);
    time_new(jj)=toc;
    %[valweights, Q1] = initialise_Hermite_transform_old(N,x);
    
    %% Some tests
    
    % Error in matrix norm

    error_matrixnorm(jj)=norm(Q'.*(1./d)' - Hforward);
    

    %norm(Q1.*valweights' - Hforward)
    
    % Error in function approximation and return
    f=@(x) 1./(2+cos(x)).*exp(-x.^2/2);
    %f=@(x) sum(exp(i*x*rand(1,3)),2).*exp(-x.^2/2);
    f_val=f(x);
    
    f_coeff1=Hforward*f_val;
    f_coeff2=Q'*((1./d).*f_val);
    
    f_val1=Hbackward*f_coeff1;
    f_val2=d.*Q*(f_coeff2);
    
    error_fn_approx_direct(jj)=norm(f_val-f_val1);
    error_fn_approx_new(jj)=norm(f_val-f_val2);

    cond_newforward(jj)=cond(Q'.*(1./d)');
    cond_Hforward(jj)=cond(Hforward);
    cond_Hbackward(jj)=cond(Hbackward);
    max_d(jj)=min(abs(d));
    min_d(jj)=max(abs(d));
    
    % To be continued
    
    
    save(strcat("data/numerical_eval_N_",num2str(Nfactor),"-",num2str(jj_max*Nfactor),".mat"))

end