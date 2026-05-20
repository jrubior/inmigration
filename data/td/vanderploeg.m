function res = vanderploeg(y,S,A,a)
% PURPOSE: Balancing by means of QL optimization (LS estimation)
% ------------------------------------------------------------------------
% SYNTAX: res = vanderploeg(y,S,A,a);
% ------------------------------------------------------------------------
% OUTPUT: res: a structure with ...
%           z       : kx1 vector of balanced variables
%           Sz      : kxk VCV of final (balanced) estimates
%           lambda  : mx1 Lagrange multipliers
% ------------------------------------------------------------------------
% INPUT:  y     : kx1 vector of unbalanced variables (initial estimates)
%         S     : kxk VCV of initial estimates
%         A     : kxm matrix of linear constraints
%         a     : 1xm vector of autonomous terms related to linear constraints
% Note: a is optional. If it is not explicitly included, the function
% assumes a=0.
% ------------------------------------------------------------------------
% LIBRARY:
% ------------------------------------------------------------------------
% SEE ALSO: bal, ras
% ------------------------------------------------------------
% REFERENCE: Van der Ploeg, F. (1982)"Reliability and the adjustment 
% of sequences of large economic accounting matrices", Journal of 
% the Royal Statistical Society, series A, vol. 145, n. 2, p. 169-194.

% written by:
%  Enrique M. Quilis

% Version 2.2 [December 2018]

% Data dimension
[k, m] = size(A);

% ------------------------------------------------------------------------
% Discrepancy
if (nargin == 3)
    a = zeros(1,m);
end
    
% Discrepancy
disc = A' * y - a' ;    
AUX = inv(A'*S*A);

% ------------------------------------------------------------------------
% Lagrange multipliers
lambda = AUX * disc;

% ------------------------------------------------------------------------
% LS balanced estimation
% ------------------------------------------------------------------------
% Levels
z = y - S * A * AUX * disc;  

% VCV
Sz = S - S * A * AUX * A' * S;

% ------------------------------------------------------------------------
%  LOADING STRUCTURE
% ------------------------------------------------------------------------
res.z      = z;
res.Sz     = Sz;
res.lambda = lambda;
