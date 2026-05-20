% PURPOSE: demo of bfl()
%          Temporal disaggregation using the Boot-Feibes-Lisman method
%---------------------------------------------------
% USAGE: bfl_d
%---------------------------------------------------
close all; clear all; clc;
% Low-frequency data: Denton's benchmark Y
load denton
% ---------------------------------------------
% Inputs
% Type of aggregation
ta = 1;   
% Minimizing the volatility of d-differenced series
d = 2;
% Frequency conversion 
sc = 12;    
% ---------------------------------------------
% Calling the function: output is loaded in a structure called res
res = bfl(Y,ta,d,sc);
% ---------------------------------------------
% Outputs
% Printed output
file_out = 'bfl.out';
tdprint(res,file_out);
edit bfl.out;
% Graphics
tdplot(res);
