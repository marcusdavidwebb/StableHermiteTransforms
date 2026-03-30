%% This is the plotting script associated with performance_evals_run.m
clear all
clc

%% Compare direct assembly of Hermite transform against stable one
addpath('modules/quadrature/')
addpath('modules/')

load("data/numerical_eval_N_multiple_transforms.mat")

%% Plot 1: Assembly time
clear h
figure(1)
h(3)=loglog(N_vec,time_GW,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(2)=loglog(N_vec,time_B,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=loglog(N_vec,time_direct,'--','Color','#c74cb9','LineWidth',3,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
h(4)=loglog(N_vec,N_vec.^2*1e-6,'-.','Color','red','LineWidth',1,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
set(gca,'FontSize',16)
ylabel('CPU time (sec)','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

legend(h,'Direct', 'Bunck', 'Golub--Welsch','$\mathcal{O}(N^2)$','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
exportgraphics(gcf,strcat("images/assembly_time_multiple.pdf"),'ContentType','vector')


%% Plot 2: Error in T
clear h
figure(2)
h(3)=loglog(N_vec,T_error_GW./T_norm_exact,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(2)=loglog(N_vec,T_error_B./T_norm_exact,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=loglog(N_vec,T_error_direct./T_norm_exact,'--','Color','#c74cb9','LineWidth',3,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
set(gca,'FontSize',16)
ylabel('$\|T-T_{\mathrm{approx}}\|_2/\|T\|_2$','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

legend(h,'Direct', 'Bunck', 'Golub--Welsch','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
ylim([1e-15,1])
exportgraphics(gcf,strcat("images/error_T_multiple.pdf"),'ContentType','vector')

%% Plot 3: Error in Tinv
clear h
figure(3)
h(3)=loglog(N_vec,Tinv_error_GW./T_inv_norm_exact,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(2)=loglog(N_vec,Tinv_error_B./T_inv_norm_exact,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=loglog(N_vec,Tinv_error_direct./T_inv_norm_exact,'--','Color','#c74cb9','LineWidth',3,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
set(gca,'FontSize',16)
ylabel('$\|T^{-1}-T_{\mathrm{approx}}^{-1}\|_2/\|T^{-1}\|_2$','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

ylim([1e-15,1])
legend(h,'Direct', 'Bunck', 'Golub--Welsch','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
exportgraphics(gcf,strcat("images/error_Tinv_multiple.pdf"),'ContentType','vector')

%% Plot 4: Condition number of T
clear h
figure(4)
h(3)=loglog(N_vec,cond_GW,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(2)=loglog(N_vec,cond_B,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=loglog(N_vec,cond_direct,'--','Color','#c74cb9','LineWidth',3,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
set(gca,'FontSize',16)
ylabel('$\mathrm{cond}(T)$','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

legend(h,'Direct', 'Bunck', 'Golub--Welsch','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
exportgraphics(gcf,strcat("images/cond_T_multiple.pdf"),'ContentType','vector')

%% Plot 5: Condition number of Tinv
clear h
figure(5)
h(3)=loglog(N_vec,cond_inv_GW,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(2)=loglog(N_vec,cond_inv_B,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=loglog(N_vec,cond_inv_direct,'--','Color','#c74cb9','LineWidth',3,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
set(gca,'FontSize',16)
ylabel('$\mathrm{cond}(T)$','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

legend(h,'Direct', 'Bunck', 'Golub--Welsch','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
ylim([1,1e16])
hold off
exportgraphics(gcf,strcat("images/cond_inv_T_multiple.pdf"),'ContentType','vector')


%% Plot 6: Accuracy of asymptotic algorithm

load("data/numerical_eval_N_multiple_transforms_asymptotic.mat")
clear h
figure(6)
h(2)=loglog(N_vec,T_error_GW./T_norm_exact,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(3)=loglog(N_vec,T_error_B./T_norm_exact,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=loglog(N_vec,N_vec.^(-4)*0.1,'--','Color','red','LineWidth',1.5,'MarkerFaceColor','white')
% loglog(N_vec,T_error_B./T_norm_exact,'-','Color','green','LineWidth',3,'MarkerFaceColor','white')
%loglog(N_vec,N_vec.^(-4),'-.','Color','red','LineWidth',1.5,'MarkerFaceColor','white')
set(gca,'FontSize',16)
ylabel('$\|T-T_{\mathrm{approx}}\|_2/\|T\|_2$','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

legend(h,'$\mathcal{O}(N^{-4})$', 'Asymptotic','Bunck','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
ylim([1e-15,1])
xlim([min(N_vec),max(N_vec)])
exportgraphics(gcf,strcat("images/error_T_asymptotic.pdf"),'ContentType','vector')


%% Plot 7: Accuracy of asymptotic algorithm in d

load("data/numerical_eval_N_multiple_transforms_asymptotic.mat")
clear h
figure(7)
h(2)=loglog(N_vec,d_error_GW,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(1)=loglog(N_vec,N_vec.^(-4)*0.1,'--','Color','red','LineWidth',1.5,'MarkerFaceColor','white')
loglog(N_vec,T_error_B./T_norm_exact,'-','Color','green','LineWidth',3,'MarkerFaceColor','white')
%loglog(N_vec,N_vec.^(-4),'-.','Color','red','LineWidth',1.5,'MarkerFaceColor','white')
set(gca,'FontSize',16)
ylabel('$\|d-d_{\mathrm{approx}}\|_2/\|d\|_2$','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

legend(h,'$\mathcal{O}(N^{-4})$', 'Asymptotic','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
ylim([1e-15,1])
xlim([min(N_vec),max(N_vec)])
exportgraphics(gcf,strcat("images/error_d_asymptotic.pdf"),'ContentType','vector')


%% Plot 8: Accuracy of asymptotic algorithm

load("data/numerical_eval_N_multiple_transforms_asymptotic.mat")
clear h
figure(8)
h(2)=loglog(N_vec,Tinv_error_GW./T_inv_norm_exact,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(3)=loglog(N_vec,Tinv_error_B./T_norm_exact,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=loglog(N_vec,N_vec.^(-4)*0.1,'--','Color','red','LineWidth',1.5,'MarkerFaceColor','white')
%loglog(N_vec,N_vec.^(-4),'-.','Color','red','LineWidth',1.5,'MarkerFaceColor','white')
set(gca,'FontSize',16)
ylabel('$\|T-T_{\mathrm{approx}}\|_2/\|T\|_2$','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

legend(h,'$\mathcal{O}(N^{-4})$', 'Asymptotic','Bunck','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
ylim([1e-15,1])
xlim([min(N_vec),max(N_vec)])
exportgraphics(gcf,strcat("images/error_T_inv_asymptotic.pdf"),'ContentType','vector')


%% Plot 9: Depiction of Q and d
figure(9)

N = 100;
[d,Q] = initialise_Hermite_transform_Golub_Welsch(N);
x = hermpts(N);

subplot(1,2,1)
m = max(abs(Q(:)));
imagesc(0:N-1,0:N-1,Q,[-m,m])
axis square
axis([0,N-1,0,N-1])
colormap turbo
colorbar
set(gca,'FontSize',16)
xlabel('$j$','Interpreter','latex', 'FontSize', 22)
ylabel('$k$','Interpreter','latex', 'FontSize', 22)
title('$Q_{kj}$','Interpreter','latex', 'FontSize', 22)

subplot(1,2,2)
plot(x,d, 'linewidth',2)
axis square
axis([1.1*x(1),1.1*x(N),0,1.1*max(d)])
grid on
set(gca,'FontSize',16)
xlabel('$x_j$','Interpreter','latex', 'FontSize', 22)
ylabel('$D_{jj}$','Interpreter','latex', 'FontSize', 22)

exportgraphics(gcf,strcat("images/depict_Q_and_d.pdf"),'ContentType','vector')