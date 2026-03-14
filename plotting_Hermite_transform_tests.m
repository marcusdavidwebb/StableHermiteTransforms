clear all
rng(1) % Fix random seed for reproducability
%% Compare direct assembly of Hermite transform against stable one

addpath('modules/quadrature/')
addpath('modules/')

% try
load("data/numerical_eval_N_10-900.mat")

%% Plot 1: Assembly time

figure(1)
semilogy(N_save, time_unscaled,'v-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
hold on
semilogy(N_save, time_new,'o-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
semilogy(N_save, time_large,'x-','linewidth',2, 'MarkerSize',10,'Color','blue','MarkerFaceColor','white')
set(gca,'FontSize',16)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)
xlabel('CPU time (sec)','Interpreter','latex', 'FontSize', 22)

legend('Unscaled', 'Scaled', 'Large N','Interpreter','latex', 'FontSize', 20,'Position',[0.27 0.8 0.1 0.1])
grid on
hold off
set(gcf, 'Position',  [100, 100, 700, 600])

exportgraphics(gcf,strcat("images/assembly_time_",num2str(Nfactor),"-",num2str(jj_max*Nfactor),".pdf"),'ContentType','vector')


%% Plot 2: error in matrix norm

figure(2)

loglog(N_save, error_matrixnorm,'<-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
hold on
loglog(N_save, error_inverse,'>-','linewidth',2, 'MarkerSize',10,'Color','blue','MarkerFaceColor','white')
loglog(N_save, 1e-3*N_save.^(-3),'--','linewidth',2, 'MarkerSize',10,'Color','red','MarkerFaceColor','white')
loglog(N_save, 1e-3*N_save.^(-4),'--','linewidth',2, 'MarkerSize',10,'Color','red','MarkerFaceColor','white')
set(gca,'FontSize',16)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)
ylabel('Transform error','Interpreter','latex', 'FontSize', 22)
legend('Forward','Inverse','Interpreter','latex', 'FontSize', 16,'Position',[0.68 0.14 0.1 0.25])
grid on
hold off
set(gcf, 'Position',  [100, 100, 700, 600])

exportgraphics(gcf,strcat("images/matrix_error_",num2str(Nfactor),"-",num2str(jj_max*Nfactor),".pdf"),'ContentType','vector')

%% Plot 3: Error in approximation
figure(3)
semilogy(N_save, error_fn_approx_unscaled,'v-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
hold on
semilogy(N_save, error_fn_approx_new,'o-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
% semilogy(N_save, error_fn_approx_unscaled,'v-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
semilogy(N_save, error_fn_approx_new,'x-','linewidth',2, 'MarkerSize',10,'Color','blue','MarkerFaceColor','white')
set(gca,'FontSize',16)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)
xlabel('$\|f-H^{-1}(H(f))\|$','Interpreter','latex', 'FontSize', 22)

legend('Unscaled', 'Scaled', 'Large N','Interpreter','latex', 'FontSize', 20,'Position',[0.27 0.8 0.1 0.1])
grid on
hold off
set(gcf, 'Position',  [100, 100, 700, 600])

exportgraphics(gcf,strcat("images/approximation_errors_",num2str(Nfactor),"-",num2str(jj_max*Nfactor),".pdf"),'ContentType','vector')

%% Figure 4
figure(4)
semilogy(N_save, cond_unscaledforward,'v-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
hold on
semilogy(N_save, cond_newforward,'o-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')


set(gca,'FontSize',16)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)
legend("cond(Q1'.*(1./d1)')","cond(Q'.*(1./d)')",'Interpreter','latex', 'FontSize', 20,'Position',[0.24 0.8 0.1 0.1])
grid on
hold off
set(gcf, 'Position',  [100, 100, 700, 600])

exportgraphics(gcf,strcat("images/condition_numbers_",num2str(Nfactor),"-",num2str(jj_max*Nfactor),".pdf"),'ContentType','vector')
% end

%% Figure 5
figure(5)
semilogy(N_save, orthogonality_unscaled,'v-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
hold on
semilogy(N_save, orthogonality_new,'o-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
semilogy(N_save, orthogonality_large,'x-','linewidth',2, 'MarkerSize',10,'Color','blue','MarkerFaceColor','white')


set(gca,'FontSize',16)
xlabel('$N$','Interpreter','latex', 'FontSize', 22)
ylabel('$\|Q^TQ-I\|$','Interpreter','latex', 'FontSize', 22)
legend("Unscaled","Scaled", "Large N",'Interpreter','latex', 'FontSize', 20,'Position',[0.24 0.8 0.1 0.1])
grid on
hold off
set(gcf, 'Position',  [100, 100, 700, 600])

exportgraphics(gcf,strcat("images/orthogonality_error_",num2str(Nfactor),"-",num2str(jj_max*Nfactor),".pdf"),'ContentType','vector')
% end