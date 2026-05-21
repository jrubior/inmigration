function [Bdraw, Sigmadraw, Rdraw_vec, Qdraw, Xdraw_vec, cholSigmadraw] = iid_sampler(PpsiTilde, cholOomegaTilde, inv_cholPphiTilde_prime, nnuTilde, setup)
% --- generate B,Sigma,Q from the posterior distribution
% Sigmadraw     = iwishrnd(PpsiTilde,nnuTilde);

% Sigma and R
Rdraw = inv_cholPphiTilde_prime * randn(setup.nvar, nnuTilde);
Sigmadraw = eye(setup.nvar) / (Rdraw * Rdraw');
Rdraw_vec = Rdraw(:);

% Bdraw
% cholSigmadraw = chol(Sigmadraw)';
% Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(setup.m*setup.nvar,1) + reshape(mmuTilde,setup.nvar*setup.m,1);
% Bdraw         = reshape(Bdraw,setup.nvar*setup.nlag+setup.nex,setup.nvar);
Bdraw = cholOomegaTilde*randn(setup.m,setup.nvar)*chol(Sigmadraw) + PpsiTilde;

% Q and X
[Qdraw,~,Xdraw]  = DrawQ(setup.nvar);
Xdraw_vec = Xdraw(:);