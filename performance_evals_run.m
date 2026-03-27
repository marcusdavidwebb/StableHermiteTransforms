%% Performance metrics of three approaches to computing the Hermite transform
clear all
rng(1) % Fix random seed for reproducability
%% Compare direct assembly of Hermite transform against stable one

addpath('modules/quadrature/')
addpath('modules/')

%% Initialise performance measures
N_vec=[4,6,8,11,16,23,32,45,64,91,128,181,256,362,512,724,1024,1448,2048,2896]; % transform sizes for experiments
zerovec=zeros(size(N_vec));

time_GW=zerovec;
time_direct=zerovec;
time_B=zerovec;

T_error_GW = zerovec;
T_error_B = zerovec;
T_error_direct = zerovec;
Tinv_error_GW = zerovec;
Tinv_error_B = zerovec;
Tinv_error_direct = zerovec;

cond_GW = zerovec;
cond_B = zerovec;
cond_direct = zerovec;

cond_inv_GW = zerovec;
cond_inv_B = zerovec;
cond_inv_direct = zerovec;

%% Evaluate the performance tests for various transform sizes
for jj=1:max(size(N_vec))
    jj
    %% Size of transform
    N=N_vec(jj);
    
    % Unstable direct algorithm
    tic
    for l=1:20 % For fairer timings repeat experiment twice
        [T, Tinv] = initialise_Hermite_transform_unstable_v2(N);
    end
    time_direct(jj)=toc/20;

    % New algorithm
    tic
    for l=1:20 % For fairer timings repeat experiment twice
        [d_GW, Q_GW] = initialise_Hermite_transform_Golub_Welsch(N);
    end
    time_GW(jj)=toc/20;
    
    % Bunck's algorithm
    tic
    for l=1:20 % For fairer timings repeat experiment twice
        [d_B, Q_B] = initialise_Hermite_transform_Bunck(N);
    end
    time_B(jj)=toc/20;
    
    % Evaluate error in approximation of T, Tinv
    % Load reference value from csv file
    d = readmatrix(strcat('data/ds_and_Qs/d_',num2str(jj),'.csv'));
    d = d(:,1);
    Q = readmatrix(strcat('data/ds_and_Qs/Q_',num2str(jj),'.csv'));

    T_error_GW(jj) = norm(diag(d_GW) * Q_GW' - diag(d) * Q');
    T_error_B(jj) = norm(diag(d_B) * Q_B' - diag(d) * Q');
    T_error_direct(jj) = norm(T - diag(d) * Q');

    Tinv_error_GW(jj) = norm(Q_GW * diag(1./d_GW) - Q * diag(1./d));
    Tinv_error_B(jj) = norm(Q_B * diag(1./d_B) - Q * diag(1./d));
    Tinv_error_direct(jj) = norm(Tinv - Q * diag(1./d));

    % Evaluate the condition number of T
    cond_GW(jj) = cond(diag(d_GW) * Q_GW');
    cond_B(jj) = cond(diag(d_B) * Q_B');
    cond_direct(jj) = cond(T);

    % Evaluate the condition number of T
    cond_inv_GW(jj) = cond(Q_GW * diag(1./d_GW));
    cond_inv_B(jj) = cond(Q_B * diag(1./d_B));
    cond_inv_direct(jj) = cond(Tinv);

    % Save the results
    save(strcat("data/numerical_eval_N_multiple_transforms.mat"),"N_vec","cond_direct","cond_B","cond_GW","cond_inv_direct","cond_inv_B","cond_inv_GW","Tinv_error_direct","Tinv_error_B","Tinv_error_GW","T_error_direct","T_error_B","T_error_GW","time_B","time_GW","time_direct")
end