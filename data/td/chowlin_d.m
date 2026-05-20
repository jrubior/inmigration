% PURPOSE: demo of chowlin()
%          Temporal disaggregation with indicators.
% 		   Chow-Lin method
%---------------------------------------------------
% USAGE: chowlin_d
%---------------------------------------------------

close all; clear all; clc;

% Loading data
load bournay_laroque;

% ---------------------------------------------
% Inputs

% Type of aggregation
ta = 1;   
% Frequency conversion 
sc = 4;    
% Method of estimation
type = 1;
% Intercept
opC = -1;
% Interval of rho for grid search
% rl = []; %Default: search on [0.05 0.99] with 100 grid points
% rl = 0.57; %Fixed value 
rl = [0.50 0.999999999 1000];
% Calling the function: output is loaded in a structure called res
res = chowlin(Y,x,ta,sc,type,opC,rl);
% Printed output
file_out = 'chowlin.out';   
tdprint(res,file_out);
edit chowlin.out;
% Graphs
tdplot(res);
