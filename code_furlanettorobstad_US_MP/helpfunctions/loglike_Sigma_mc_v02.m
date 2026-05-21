function loglik_prop = loglike_Sigma_mc_v02(z_prop,function_restrictions,Bdraw,Qdraw,chol_lower_OomegaTilde,chol_lower_OomegaTilde_inv,setup)
% likelihood eval y given z_prop
% z_prop = [ddelta, ggamma]


% ---

R_prop = reshape(z_prop,setup.nvar,setup.nnuTilde);


% Sigma_prop=inv(R_prop*R_prop');
Sigma_prop=eye(setup.nvar)/(R_prop*R_prop');


S_prop = function_restrictions(Bdraw, Sigma_prop, Qdraw,setup);


chol_lower_Sigma_prop = chol(Sigma_prop,'lower');
% chol_lower_OomegaTilde;


% ---
% p(B|Sigma)

% ---
% Fast likelihood calculation using kronecker product structure
% note that
% chol_upper(V) = kron(chol_upper(Sigma),chol_upper(Omega))
% Want to find xRinv such that  xRinv * chol_upper(V) = X0 where xRinv
% and X0 are row vectors
% Then, chol_lower(V) * xRinv' = X0'
% Then, chol_lower(Omega) * unvec(xRinv') * chol_lower(Sigma)' = unvec(X0')
% And therefore, unvec(xRinv') = chol_lower(Omega) \ unvec(X0') / chol_lower(Sigma)'
% dim(Omega) = n_Omega (setup.m)
% dim(Sigma) = n_Sigma (setup.nvar)
% ---

% % % actual calculation
% unvec_X0_prime = reshape( ( vec(Bdraw) - vec(setup.PpsiTilde) )', setup.m, setup.nvar);
% % unvec_xRinv_prime = chol_lower_OomegaTilde \ unvec_X0_prime / chol_lower_Sigma_prop;
% unvec_xRinv_prime = chol_lower_OomegaTilde_inv * unvec_X0_prime / chol_lower_Sigma_prop;
% xRinv_prime = unvec_xRinv_prime(:);
% 
% logSqrtDetSigma = sum(log(diag(chol_lower_Sigma_prop))) * setup.nvar + sum(log(diag(chol_lower_OomegaTilde))) * setup.m;
% 
% quadform = sum(xRinv_prime.^2, 1); %xRinv_prime is column vector
% 
% logpdfBgivenSigma = (-0.5*quadform - logSqrtDetSigma - size(xRinv_prime,1)*log(2*pi)/2);

% % actual calculation: jonas
unvec_X0_prime = reshape( ( vec(Bdraw) - vec(setup.PpsiTilde) ), setup.m, setup.nvar);
% unvec_xRinv_prime = chol_lower_OomegaTilde \ unvec_X0_prime / chol_lower_Sigma_prop;
%unvec_xRinv_prime = chol_lower_OomegaTilde_inv' * unvec_X0_prime / chol_lower_Sigma_prop;% minchul found typo
unvec_xRinv_prime = chol_lower_OomegaTilde_inv * unvec_X0_prime / chol_lower_Sigma_prop';
xRinv_prime = unvec_xRinv_prime(:);

logSqrtDetSigma = sum(log(diag(chol_lower_Sigma_prop))) * setup.m + sum(log(diag(chol_lower_OomegaTilde))) * setup.nvar;

quadform = sum(xRinv_prime.^2, 1); %xRinv_prime is column vector

logpdfBgivenSigma = (-0.5*quadform - 1*logSqrtDetSigma - size(xRinv_prime,1)*log(2*pi)/2);




loglik_prop = -Inf;
if S_prop==1
    
loglik_prop = logpdfBgivenSigma ;
end
