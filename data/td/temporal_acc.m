function [za] = temporal_acc(z,op1,sc)
% PURPOSE: Annual cumulative to date (YTD) high frequency time series
% ------------------------------------------------------------
% SYNTAX: za = temporal_acc(z,op1,sc);
% ------------------------------------------------------------
% OUTPUT: za: nxk YTD accumulated vector time series
% ------------------------------------------------------------
% INPUT:  z: nxk ---> vector of high frequency data
%         op1: type of temporal aggregation 
%         op1=1 ---> sum (flow)
%         op1=2 ---> average (index)
%         op1=3 ---> last element (stock) ---> interpolation
%         op1=4 ---> first element (stock) ---> interpolation
%         sc: number of high frequency data points 
%            for each low frequency data points
% ------------------------------------------------------------
% LIBRARY: acc
% ------------------------------------------------------------
% SEE ALSO: temporal_agg, copylow, ssampler, aggreg_p
% ------------------------------------------------------------

% written by:
%  Enrique M. Quilis

% Version 2.1 [December 2018]

% Data dimension
n = length(z);

% Computing implicit numer of low-frequency data
N = fix(n/sc);

% Generating accumulation matrix
A = acc(op1,n,sc);

% Computing accumulated time series
za = A * z;

