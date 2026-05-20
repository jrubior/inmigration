function res = uni_chowlin(Y,ta,sc,opx)
% PURPOSE: Temporal disaggregation using the Chow-Lin method without
% indicator (i.e., passing time t=1..n as input to chowlin()).
% ------------------------------------------------------------
% SYNTAX: res = uni_chowlin(Y,ta,sc);
% ------------------------------------------------------------
% OUTPUT: res: a structure
%           res.meth    ='Chow-Lin';
%           res.ta      = type of disaggregation
%           res.type    = method of estimation
%           res.opC     = option related to intercept
%           res.N       = nobs. of low frequency data
%           res.n       = nobs. of high-frequency data
%           res.pred    = number of extrapolations
%           res.sc      = frequency conversion between low and high freq.
%           res.p       = number of regressors (including intercept)
%           res.Y       = low frequency data
%           res.x       = high frequency indicators
%           res.y       = high frequency estimate
%           res.y_dt    = high frequency estimate: standard deviation
%           res.y_lo    = high frequency estimate: sd - sigma
%           res.y_up    = high frequency estimate: sd + sigma
%           res.u       = high frequency residuals
%           res.U       = low frequency residuals
%           res.beta    = estimated model parameters
%           res.beta_sd = estimated model parameters: standard deviation
%           res.beta_t  = estimated model parameters: t ratios
%           res.rho     = innovational parameter
%           res.sigma_a = variance of shocks
%           res.aic     = Information criterion: AIC
%           res.bic     = Information criterion: BIC
%           res.val     = Objective function used by the estimation method
%           res.wls     = Weighted least squares as a function of rho
%           res.loglik  = Log likelihood as a function of rho
%           res.r       = grid of innovational parameters used by the estimation method
% ------------------------------------------------------------
% INPUT: Y: Nx1 ---> vector of low frequency data
%        ta: 1x1 -> type of disaggregation
%            ta=1 ---> sum (flow)
%            ta=2 ---> average (index)
%            ta=3 ---> last element (stock) ---> interpolation
%            ta=4 ---> first element (stock) ---> interpolation
%        sc: 1x1 -> number of high frequency data points for each low frequency data points 
%            Some examples:
%            sc= 4 ---> annual to quarterly
%            sc=12 ---> annual to monthly
%            sc= 3 ---> quarterly to monthly
%        opx: 1x1 -> model for the indicator: 1 => linear trend, 
%                    2 => quadratic trend [default]
% ------------------------------------------------------------
% LIBRARY: chowlin
% ------------------------------------------------------------
% SEE ALSO: bfl
% ------------------------------------------------------------
% REFERENCE: Chow, G. and Lin, A.L. (1971) "Best linear unbiased 
% distribution and extrapolation of economic time series by related 
% series", Review of Economic and Statistics, vol. 53, n. 4, p. 372-375.
% Bournay, J. and Laroque, G. (1979) "Reflexions sur la methode 
% d'elaboration des comptes trimestriels", Annales de l'INSEE, n. 36, p. 3-30.

% ------------------------------------------------------------
% written by:
%  Enrique M. Quilis

% Version 1.1 [December 2018]

% Determining default options (opx=2 => quadratic trend)
if (nargin == 3)
    opx = 2; 
end

% Data dimension.
N = length(Y);
n = sc*N;

% Generating time as "explanatory" high-frequency indicator.
switch opx
    case 1 %Linear trend
        x = (1:n)';
    case 2 %Quadratic trend
        x1 = (1:n)';
        x2 = x1 .^2;
        x = [x1 x2];
        clear x1 x2;
end

% Method of estimation: MLE.
type = 1;

% Intercept: pretest.
opC = -1;

% Interval of rho for grid search: default.
rl = [];

% Calling the Chow-Lin function: output is loaded in a structure res.
res = chowlin(Y,x,ta,sc,type,opC,rl);

