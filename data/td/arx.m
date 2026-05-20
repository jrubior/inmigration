function res = arx(rex)
% PURPOSE: Unrestricted estimation of an univariate ARX model
% ------------------------------------------------------------
% SYNTAX:  res = arx(rex)
% ------------------------------------------------------------
% OUTPUT: a structure ...
%         res.meth       = 'ARX MODEL: UNRESTRICTED ESTIMATION';
%         res.parm       = [n,nobs,k,p];
%         res.mu         = [mu]; <========= Intercept
%         res.mu_dt      = [mu_dt];
%         res.mu_t       = [mu_t];
%         res.b          = [b];  <========= Exogenous variables parameters
%         res.b_dt       = [b_dt];
%         res.b_t        = [b_t];
%         res.phi        = [phi]; <========= Lagged variables parameters
%         res.phi_dt     = [phi_dt];
%         res.phi_t      = [phi_t];
%         res.u          = [u];   <========= Residuals
%         res.sigma      = [sigma];
%         res.perform    = [log_like,aic,bic];
%         res.G              = VCV matrix of estimated parameters
% ------------------------------------------------------------
% INPUT: a structure ...
%        X    : nxm matrix of exogenous series, columnwise
%        z    : nx1 vector of original series, columnwise
%        opC  : 1x1 indicator to include (1) an intercept or not (0)
%        p    : 1x1 order of the model
% ------------------------------------------------------------
% LIBRARY: lag, uisolator
% ------------------------------------------------------------
% SEE ALSO: barx, arx_print

% ------------------------------------------------------------
% written by:
%  Enrique M. Quilis

% January 2009
% Version 1.2

% ------------------------------------------------------------
% Loading the structure

% Exogenous variables
Xx = rex.X;

% Intercept: yes (1) or not (0)
opC = rex.opC;

% Basic variable
z = rex.z;

% Order of AR model
p = rex.p;

% ------------------------------------------------------------
% Computing basic parameters

[n,k]   = size(z);
if (k ~= 1)
    error('*** ONLY UNIVARIATE MODELING, PLEASE ***');
end
    
[aux,m] = size(Xx);

nobs = n-p;  % Number of effective observations

% ------------------------------------------------------------
%       PREPARING THE MODEL IN LINEAR MODEL FORMAT: Y=XB+U
% ------------------------------------------------------------

Xz = lag(z,p);
Xz = Xz(:,2:end);  % Dropping vector of ones (no intercept is assumed)

z  = z(p+1:end);

if (isempty(Xx) == 0)  % Exogenous variables are included
   Xx = Xx(p+1:end,:);
end

% Intercept

switch opC
case 0 % No intercept
   X = [Xx Xz];
case 1 % Intercept
   X = [ones(nobs,1) Xx Xz];
end

% ------------------------------------------------------------
%       ESTIMATION OF THE MODEL: OLS

beta = (X'*X) \ (X'*z);

% Residuals   

u = z - X * beta ;
 
% variance of residuals   

sigma = var(u,1);

% VCV and t-ratios of beta parameters
G           = sigma * inv(X'*X);
beta_dt  = sqrt(diag(G));
beta_t    = beta ./ beta_dt;

% Performance measures

log_like = -0.5*nobs*(1*(1+log(2*pi))+log(sigma));

npar = opC + m + p;
aic = -(2*log_like/nobs) + (2*npar)/nobs;
bic = -(2*log_like/nobs) + npar*(log(nobs)/nobs);

% ------------------------------------------------------------
%       ISOLATING THE MATRICES
% ------------------------------------------------------------

[mu,b,phi] = uisolator(beta,opC,m,p);
[mu_dt,b_dt,phi_dt] = uisolator(beta_dt,opC,m,p);
[mu_t,b_t,phi_t] = uisolator(beta_t,opC,m,p);

% ------------------------------------------------------------
%       LOADING THE STRUCTURE

res.meth       = 'ARX MODEL: UNRESTRICTED ESTIMATION';
res.parm       = [n,nobs,k,p];
res.mu         = [mu];
res.mu_dt      = [mu_dt];
res.mu_t       = [mu_t];
res.b          = [b];
res.b_dt       = [b_dt];
res.b_t        = [b_t];
res.phi        = [phi];
res.phi_dt     = [phi_dt];
res.phi_t      = [phi_t];
res.u          = [u];
res.sigma      = [sigma];
res.perform    = [log_like,aic,bic];
res.G             = [G];
   
   
