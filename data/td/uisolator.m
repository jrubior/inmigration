function [mu,b,phi] = uisolator(beta,opC,m,p)
% PURPOSE: Creates a sequence of vectors of dimension 1 (mu)
%          m (beta for exog.) and p (lagged endogenous)
% ------------------------------------------------------------
% SYNTAX: [mu,b,phi] = uisolator(beta,opC,m,p);
% ------------------------------------------------------------
% OUTPUT: Intercept mu : 1x1, b : mx1, phi : px1
%         Note thar mu and/or b may be empty 
% ------------------------------------------------------------
% INPUT: beta: vector to be segmented
%        opC : 1x1 indicator to include (1) an intercept or not (0)
%        m   : number of exogenous variables
%        p   : number of lags
% ------------------------------------------------------------
% LIBRARY: 
% ------------------------------------------------------------
% SEE ALSO: isolator 

% ------------------------------------------------------------
% written by:
%  Enrique M. Quilis

% September 2007
% Version 1.1

% ------------------------------------------------------------
% Isolating intercept

switch opC
case 0 % No intercept
   mu = [];
case 1
   mu = beta(1,1);
end

% ------------------------------------------------------------
% Isolating beta for exogenous variables

switch m
case 0 % No exogenous variables
   b = [];
otherwise % Exogenous variables
   b = beta(2:m+1,1);
end

% ------------------------------------------------------------
% Isolating beta for lagged endogenous variable

phi = beta(m+opC+1:end,1);


