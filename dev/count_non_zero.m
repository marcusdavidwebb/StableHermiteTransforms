clear all
rng(1) % Fix random seed for reproducability
%% Add paths

addpath('modules/quadrature/')
addpath('modules/')

%% Initialise
N_vec=[4,6,8,11,16,23,32,45,64,91,128,181,256,362,512,724,1024,1448,2048,2896,4096]; % transform sizes for experiments
zerovec=zeros(size(N_vec));

count=zerovec;

%% Evaluate
for jj=1:max(size(N_vec))
    jj
    [d_GW, Q_GW] = initialise_Hermite_transform_Golub_Welsch(N_vec(jj));
    tol = 1e-10;              % choose a tolerance
    count(jj) = sum(abs(Q_GW(:)) < tol);
    % count(jj)=nnz(Q_GW);
end

%% Plotting
figure(1)
loglog(N_vec,count,'b-','LineWidth',2)
hold on
loglog(N_vec,N_vec.^1.5,'r-.','LineWidth',2)
loglog(N_vec,N_vec.^2,'r-.','LineWidth',2)

legend('\# $|Q_{ij}|<10^{-10}$','$N^{1.5}, N^2$','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off

exportgraphics(gcf,strcat("../dev/num_nonzeroQ.pdf"),'ContentType','vector')


%% Plotting
figure(2)
loglog(N_vec,N_vec.^2-count,'b-','LineWidth',2)
hold on
loglog(N_vec,N_vec.^1.5,'r-.','LineWidth',2)
loglog(N_vec,N_vec.^2,'r-.','LineWidth',2)

legend('\# $|Q_{ij}|>10^{-10}$','$N^{1.5}, N^2$','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off

exportgraphics(gcf,strcat("../dev/num_zeroQ.pdf"),'ContentType','vector')