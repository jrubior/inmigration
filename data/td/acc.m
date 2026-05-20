function [A] = acc(op1,n,sc)
% PURPOSE: Generate an accumulation matrix
% ------------------------------------------------------------
% SYNTAX: [A] = acc(op1,n,sc); 
% ------------------------------------------------------------
% OUTPUT: A: nxn accumulation matrix
% ------------------------------------------------------------
% INPUT:  op1: type of temporal aggregation 
%         op1=1 ---> sum (flow)
%         op1=2 ---> average (index)
%         n: number of high frequency data
%         sc: number of high frequency data points 
%            for each low frequency data points (freq. conversion)
% ------------------------------------------------------------
% SEE ALSO: temporal_agg, ssampler, copylow
% ------------------------------------------------------------

% written by:
%  Enrique M. Quilis

% Version 2.0 [December 2018]

% Computing implicit numer of low-frequency data
N = fix(n/sc);

% Counting number of observations of last year
npred = n - sc*N;

% Expanding the number of years to include open last year (if required)
if (npred > 0)
    N = N + 1;
end

% Auxiliary matrix
A1 = tril(ones(sc,sc));

% Selecting accumulation or averaging
switch op1
    case 1
        % Do nothing
    case 2
        A1 = A1 / sc;
    otherwise
        error ('*** Improper op1 ***');
end

% Accumulation matix
A = kron(eye(N),A1);

% Adjusting the accumulation matrix to the real dimension of input data
A = A(1:n,1:n);
