% PURPOSE: demo of backtest()
%          Exploring past performance
%---------------------------------------------------
% USAGE: backtest_d
%---------------------------------------------------

close all; clear all; clc;

% Loading data
load bournay_laroque;
  
% ---------------------------------------------
% Inputs
% ---------------------------------------------

% Type of aggregation
ta = 1;   
% Frequency conversion 
sc = 4;    
% Method of estimation
type = 1;
% Intercept
opC = 1;
% Interval of rho for grid search
rl = [0 0.9999999999  100];
% Temporal disaggregation
rex = chowlin(Y,x,ta,sc,type,opC,rl);
% Number of low-frequency periods in which the history of 
% revisions will be perfomed
Nh = 3;
res = backtest(rex,Nh);
% Temporal aggregates: recursive estimates
Yh = temporal_agg(res.yh,ta,sc);
% Graph
plot(res.yh)
    title('Revisions')
