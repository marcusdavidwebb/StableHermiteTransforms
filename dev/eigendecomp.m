clear all;
clc;
%% Parameters
jmax=10;
cost_sparse=zeros(jmax,1);
cost_dense=zeros(jmax,1);


%% Compare cost

for j=1:10
    n=2^j;
    a=sqrt((1:(n-1))'/2)
    A=spdiags(a,-1,n,n);
    A=A+A';
    
    [V,D]=eig(full(A));
    d=diag(D);
    tic
    for l=1:n
        [V,D]=eigs(A-d(l)*speye(n),1);
    end
    cost_sparse(j)=toc;

    tic
    for l=1:10
        [V,D]=eig(full(A));
    end
    cost_dense(j)=toc/10;


end

%% Plot results

semilogy(2.^(1:jmax),cost_sparse)
hold on
semilogy(2.^(1:jmax), cost_dense);
xlabel('Matrix Size (n)');
ylabel('Computation Time (seconds)');
legend('Sparse', 'Dense');
title('Comparison of Sparse and Dense Eigenvalue Computation');
grid on;
hold off;