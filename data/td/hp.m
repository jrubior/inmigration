function [trend, cycle]=hp(Z,lambda)
% PURPOSE: Hodrick-Prescott filtering. Two-sided. Matrix computations
% ------------------------------------------------------------
% SYNTAX: [trend, cycle] = hp(Z,lambda);
% ------------------------------------------------------------
% OUTPUT: trend: nxk --> long-term trend
%         cycle: nxk --> cycle as deviation from long-term trend
% ------------------------------------------------------------
% INPUT: Z: nxk      --> vector of observed time series
%        lambda: 1x1 --> weight in the penalty function
% ------------------------------------------------------------
% LIBRARY: dif
% ------------------------------------------------------------

% ------------------------------------------------------------
% written by:
%  Enrique M. Quilis

% Version 3.0 [October 2015]

% Size of the input data
n = size(Z,1);

% Second difference matrix
D = dif(2,n);
D(1:2,:) = [];

% Matrix form of the filter (non inverted)
H = eye(n) + lambda*(D'*D);

% Computing long-term trend
trend = H \ Z;

% Computing cycle as the deviation from observed trend
cycle = Z - trend;

% Filter in matrix form
Hhp = inv(H);
