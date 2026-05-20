function [z_acc] = accYTD(z,op1,sc)
% PURPOSE: Annual cumulative to date (YTD) high frequency time series
% ---------------------------------------------------------------------
% SYNTAX: [z_acc] = accYTD(z,op1,sc);
% ---------------------------------------------------------------------
% OUTPUT: z_acc: nxk YTD accumulated time series
% ---------------------------------------------------------------------
% INPUT:  z: nxk ---> vector of high frequency data
%         op1: type of temporal aggregation 
%         op1=1 ---> sum (flow)
%         op1=2 ---> average (index)
%         op1=3 ---> last element (stock) ---> interpolation
%         op1=4 ---> first element (stock) ---> interpolation
%         sc: number of high frequency data points 
%            for each low frequency data points
% ---------------------------------------------------------------------
% LIBRARY: copylow, temporal_agg, ssampler, aggreg_p
% ---------------------------------------------------------------------

% written by:
%  Enrique M. Quilis

% Version 1.1 [December 2018]

% ---------------------------------------------------------------------
% Data dimension
n = size(z,2);

% ---------------------------------------------------------------------
% Generating accumulation matrix
AC = acc(op1,n,sc);

% ---------------------------------------------------------------------
% Generating annually accumulated time series
z_acc = AC * z;
