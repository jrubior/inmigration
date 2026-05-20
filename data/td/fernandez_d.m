% PURPOSE: demo of fernandez()
%          Temporal disaggregation with indicators.
% 			  Fernandez method
%---------------------------------------------------
% USAGE: fernandez_d
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
% Intercept
opC = -1;
% Calling the function: output is loaded in a structure called res
res = fernandez(Y,x,ta,sc,opC);
% Printed output
file_out = 'fernandez.out';   
tdprint(res,file_out);
edit fernandez.out;
% Calling graph function
tdplot(res);
