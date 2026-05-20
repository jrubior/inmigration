% PURPOSE: demo of litterman()
%          Temporal disaggregation with indicators.
% 			  Litterman method
%---------------------------------------------------
% USAGE: litterman_d
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
rl = [-0.99 0.99 500];
% Calling the function: output is loaded in a structure called res
res = litterman(Y,x,ta,sc,type,opC,rl);
% Printed output
file_out = 'litterman.out';   
tdprint(res,file_out);
edit litterman.out;
% Graphs
tdplot(res);
