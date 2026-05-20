function [y] = temporal_agg_p(z,op1,sc)
% PURPOSE: Temporal aggregation of a time series preserving its dimension
% ------------------------------------------------------------
% SYNTAX: y = temporal_agg_p(z,op1,sc);
% ------------------------------------------------------------
% OUTPUT: y: nx1 temporally aggregated series, missing=0
% ------------------------------------------------------------
% INPUT:  z: nx1 ---> vector of high frequency data
%         op1: type of temporal aggregation 
%         op1=1 ---> sum (flow)
%         op1=2 ---> average (index)
%         op1=3 ---> last element (stock) ---> interpolation
%         op1=4 ---> first element (stock) ---> interpolation
%         sc: number of high frequency data points 
%            for each low frequency data points
%         Note: n = sc x N
% ------------------------------------------------------------
% LIBRARY: aggreg
% ------------------------------------------------------------

% written by:
%  Enrique M. Quilis

% Version 1.1 [December 2018]

% Data dimension
n = size(z,1);

% ------------------------------------------------------------
% Computes the number of low frequency points. 
% Low frequency periods should be complete
N = fix(n/sc);
C = aggreg_p(op1,n,sc);

% -----------------------------------------------------------
% Expanding the aggregation matrix to perform
% extrapolation if needed.
if (n > sc * N)
   pred = n - sc*N;           % Number of required extrapolations 
   C=[C zeros(N,pred)];
else
   pred = 0;
end

y = C*z;
