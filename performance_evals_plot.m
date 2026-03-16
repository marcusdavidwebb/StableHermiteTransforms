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
h(3)=semilogy(N_vec,time_GW,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(2)=semilogy(N_vec,time_B,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=semilogy(N_vec,time_direct,'--','Color','#c74cb9','LineWidth',3,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
set(gca,'FontSize',16)
ylabel('CPU time (sec)','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

legend(h,'Direct', 'Bunck', 'Golub--Welsch','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
exportgraphics(gcf,strcat("images/assembly_time_multiple.pdf"),'ContentType','vector')


%% Plot 2: Error in T
clear h
figure(2)
h(3)=semilogy(N_vec,T_error_GW,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(2)=semilogy(N_vec,T_error_B,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=semilogy(N_vec,T_error_direct,'--','Color','#c74cb9','LineWidth',3,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
set(gca,'FontSize',16)
ylabel('$\|T-T_{\mathrm{approx}}\|_2$','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

legend(h,'Direct', 'Bunck', 'Golub--Welsch','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
ylim([1e-15,1])
exportgraphics(gcf,strcat("images/error_T_multiple.pdf"),'ContentType','vector')

%% Plot 3: Error in Tinv
clear h
figure(3)
h(3)=semilogy(N_vec,Tinv_error_GW,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(2)=semilogy(N_vec,Tinv_error_B,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=semilogy(N_vec,Tinv_error_direct,'--','Color','#c74cb9','LineWidth',3,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
set(gca,'FontSize',16)
ylabel('$\|T^{-1}-T_{\mathrm{approx}}^{-1}\|_2$','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

ylim([1e-15,1])
legend(h,'Direct', 'Bunck', 'Golub--Welsch','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
exportgraphics(gcf,strcat("images/error_Tinv_multiple.pdf"),'ContentType','vector')

%% Plot 4: Condition number of T
clear h
figure(4)
h(3)=semilogy(N_vec,cond_GW,'-','Color','blue','LineWidth',3,'MarkerFaceColor','white')
hold on
h(2)=semilogy(N_vec,cond_B,'-.','Color','#edb120','LineWidth',3,'MarkerFaceColor','white')
h(1)=semilogy(N_vec,cond_direct,'--','Color','#c74cb9','LineWidth',3,'MarkerFaceColor','white')%,'Color','#edb120','LineWidth',2)
set(gca,'FontSize',16)
ylabel('$\mathrm{cond}(T)$','Interpreter','latex', 'FontSize', 22)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)

legend(h,'Direct', 'Bunck', 'Golub--Welsch','Interpreter','latex', 'FontSize', 16,'Location','northwest')
grid on
hold off
exportgraphics(gcf,strcat("images/cond_T_multiple.pdf"),'ContentType','vector')