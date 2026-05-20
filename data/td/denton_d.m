% PURPOSE: Demo of denton()
%          Temporal disaggregation with indicator.
% 		   Denton method, additive and proportional variants.
%---------------------------------------------------
% USAGE: denton_d
%---------------------------------------------------

close all; clear all; clc;

% Loading data
load denton;

% ---------------------------------------------
% Inputs.
% ---------------------------------------------

% Type of aggregation.
ta = 1;   
% Minimizing the volatility of d-differenced series.
d = 2;
% Frequency conversion.
sc = 4;
% Variant: 1=additive, 2=proportional.
op1 = 2;
% Calling the function: output is loaded in a structure called res.
res = denton(Y,x,ta,d,sc,op1);
% Calling printing function.
% Name of ASCII file for output.
file_out = 'denton.out';   
tdprint(res,file_out);
edit denton.out;
% Final plots
figure;
t = 1:length(res.y);
subplot(2,1,1)
plot(t,copylow(Y,2,sc),t,res.y,'r-');
    title('Low frequency input');
    legend('Benchmark','Interpolation');
    title('Denton procedure');
subplot(2,1,2)
plot(t,copylow(res.U,2,sc),t,res.u,'r-');
    legend('Low freq.','High freq.');
    title('Residuals');
