function res = ols_regress(Y,X)
% PURPOSE: OLS estimates of linear regression model Y=Xb+U
% ------------------------------------------------------------
% SYNTAX: res = ols_regress(Y,X);
% ------------------------------------------------------------
% OUTPUT: res: a structure
%         res.beta     = OLS Estimates of b: kx1 
%         res.beta_se  = Standard errors of OLS Estimates of b: kx1 
%         res.beta_t   = t-ratios of OLS Estimates of b: kx1 
% ------------------------------------------------------------
% INPUT: Y: Nx1 ---> response data
%        X: nxp ---> regressors
% ------------------------------------------------------------
% NOTE: No automatic estimation of intercept is performed
% ------------------------------------------------------------

% written by:
%  Enrique M. Quilis

% Version 1.0 [December 2018]

% OLS Estimates
beta = (X'*X) \ (X'*Y);
% Residuals
e = Y - X*beta;
% Estimate of sigma
sigma = (e'*e) / (length(Y)-size(X,2));
% VCV matrix of beta
S_beta = sigma*inv(X'*X);
% Selecting s.e. of beta
beta_se = sqrt(diag(S_beta));
% t-ratios of beta
beta_t = beta ./ beta_se;

% -----------------------------------------------------------
% Loading the structure
% -----------------------------------------------------------
res.beta = beta;
res.beta_se = beta_se;
res.beta_t = beta_t;
