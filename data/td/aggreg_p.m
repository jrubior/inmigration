function [C] = aggreg_p(op1,n,sc);
% PURPOSE: Temporal aggregation matrix preserving input dimension
% ------------------------------------------------------------
% SYNTAX: [C] = aggreg_p(op1,n,sc); 
% ------------------------------------------------------------
% OUTPUT: C: nxn temporal aggregation matrix
% ------------------------------------------------------------
% INPUT:  op1: type of temporal aggregation 
%         op1=1 ---> sum (flow)
%         op1=2 ---> average (index)
%         op1=3 ---> last element (stock) ---> interpolation
%         op1=4 ---> first element (stock) ---> interpolation
%         n: number of high frequency data
%         sc: number of high frequency data points 
%            for each low frequency data points (freq. conversion)
% ------------------------------------------------------------
% LIBRARY: aggreg_v
% ------------------------------------------------------------
% SEE ALSO: aggreg, temporal_agg
% ------------------------------------------------------------
% NOTE: Use aggreg_v_X for extended interpolation

% written by:
%  Enrique M. Quilis

% Version 2.0 [December 2018]

% ------------------------------------------------------------
% Computing implicit numer of low-frequency data
N = fix(n/sc);

% Counting number of observations of last year
npred = n - sc*N;

% Expanding the number of years to include open last year (if required)
if (npred > 0)
    N = N + 1;
end

% ------------------------------------------------------------
% Generation of temporal aggregation matrix
c = aggreg_v(op1,sc);
C1 = [zeros(sc-1,sc) ; c];
C = kron(eye(N),C1);
% Adjusting the accumulation matrix to the real dimension of input data
C = C(1:n,1:n);


