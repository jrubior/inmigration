% PURPOSE: demo of ssc()
%          Temporal disaggregation with indicators.
%          Santos-Silva & Cardoso method
%---------------------------------------------------
% USAGE: ssc_d
%---------------------------------------------------

close all; clear all; clc;

% Loading data
load ssc;

% ---------------------------------------------
% Inputs

% Type of aggregation
ta=1;
% Frequency conversion
s=4;
% Method of estimation
type=0;
% Intercept
opC = 1;
% Interval of rho for grid search
% rl = [];
% rl = 0.57;
rl = [-0.99 0.99 500];

% Note: the grid search applied in the undelying ssc procedure generates 
% a warning when phi=0. This warning is muted.
warning off MATLAB:nearlySingularMatrix
% Calling the function: output is loaded in a structure called res
res = ssc(Y,x,ta,s,type,opC,rl);
warning on MATLAB:nearlySingularMatrix

% Printed output
file_out = 'ssc.out';
tdprint(res,file_out);
edit ssc.out;

% Graphs
tdplot(res);
