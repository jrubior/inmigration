% PURPOSE: Demo of chowlin_uni()
%          Temporal disaggregation without indicator and Chow-Lin.
%---------------------------------------------------
% USAGE: chowlin_uni_d
%---------------------------------------------------

close all; clear all; clc;

% Low-frequency data: Denton's benchmark Y
load denton;
    
% ---------------------------------------------
% Inputs 
% ---------------------------------------------
% Type of aggregation.
ta = 2;   
% Frequency conversion.
sc = 3;
% Model for the indicator: 1 => linear trend, 2 => quadratic trend
opx = 2;
% Calling the function: output is loaded in a structure called res.
res = uni_chowlin(Y,ta,sc,opx);
% Calling printing function.
file_out = 'cl_uni.out';   
tdprint(res,file_out);
edit cl_uni.out;
% Graphs
t = 1:res.n;
plot(t,copylow(Y,1,sc),t,res.y,'r-',t,res.y_lo,'r:',t,res.y_up,'r:');
    legend('Benchmark','Interpolation','Interpolation: s.e.','Interpolation: s.e.','Location','best');
    title('Chow-Lin interpolation');
