function [ttheta_old,ssigma_old] = BSigmaTOtthetassigma2(Bdraw,Sigmadraw)

% mapping Bhat,SigmaHat to theta_hat = [bbeta_1,...,beta_n,alpha],sigma_hat^2
% the columns of Bhat corresponds to beta_1,...,beta_n



nvar=size(Sigmadraw,1);
ttheta_old = cell([nvar,1]);

% Standard Cholesky decomposition
C = chol(Sigmadraw, 'lower');

% Extract diagonal elements
d = diag(C);

% Create diagonal matrix D
D = diag(d.^2);

% Normalize to get L with ones on diagonal
Ainv = C ./ d';

A = Ainv\eye(size(Ainv,1));
Bdraw= Bdraw*A';

for i=1:nvar
ttheta_old{i} = Bdraw(:,i);
end


A_id = nonzeros(tril(reshape(1:nvar^2,nvar,nvar),-1)');
aalpha_all = A(A_id);
    count_alp = 0;
for i=1:nvar

    aalpha_oldi=aalpha_all(count_alp+1:count_alp+i-1);
ttheta_old{i} = [ttheta_old{i};aalpha_oldi];
        count_alp = count_alp + i-1; 
end
ssigma_old = diag(D);