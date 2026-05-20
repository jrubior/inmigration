% PURPOSE: demo of guerrero()
%          Temporal disaggregation with indicators.
%          Guerrero ARIMA-based method
%---------------------------------------------------
% USAGE: guerrero_d
%---------------------------------------------------

clear all; clc; close all;

% ------------------------------------------------------------
% Loading data
load guerrero;
    
% ---------------------------------------------
% Inputs for td library

% Type of aggregation
ta = 1;   
% Frequency conversion 
sc = 12;    
% Intercept
opC = 1;

% Model for  w: (0,1,1)(1,0,1)
rexw.ar_reg = [1];
rexw.d  = 1;
rexw.ma_reg = [1 -0.40];

rexw.ar_sea = [1 0 0 0 0 0 0 0 0 0 0 0 -0.85];
rexw.bd = 0;
rexw.ma_sea = [1 0 0 0 0 0 0 0 0 0 0 0 -0.79];

rexw.sigma = 4968.716^2;

% Model for the discrepancy: (1,2,0)(1,0,0)
% See: Martinez and Guerrero, 1995, Test, 4(2), 359-76.

rexd.ar_reg = [1 -0.43]; 
rexd.d  = 2;
rexd.ma_reg = [1];

rexd.ar_sea = [1 0 0 0 0 0 0 0 0 0 0 0 0.62]; 
rexd.bd = 0;
rexd.ma_sea = [1];

rexd.sigma = 76.95^2;

% Calling the function: output is loaded in a structure called res
res = guerrero(Y,x,ta,sc,rexw,rexd,opC);

% Printed output
file_out = 'guerrero.out';   
td_print(res,file_out);  
edit guerrero.out;

% Graphs
td_plot(res);
