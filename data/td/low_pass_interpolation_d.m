% PURPOSE: demo of low_pass_interpolation()
%          Temporal disaggregation using low-pass filtering
%          Low-pass filter = Hodrick-Prescott
% ---------------------------------------------------------------------
% USAGE: low_pass_interpolation_d
% ---------------------------------------------------------------------

clear all; close all; clc;
% ---------------------------------------------------------------------
% Annual GDP. Spain. 1960-2020. (AMECO Database).
load Spain_GDP;
Z = Y;
% ---------------------------------------------------------------------
% Sample conversion
sc = 4;
% ---------------------------------------------------------------------
% Hodrick-Prescott parameter
lambda = 1600;
% ---------------------------------------------------------------------
% Denton parameters:
% Type of aggregation
ta = 2;   
% Minimizing the volatility of d-differenced series
d = 1;
% ---------------------------------------------------------------------
% Calling function
[z,w,x] = low_pass_interpolation(Z,ta,d,sc,lambda);
% ---------------------------------------------------------------------
% Plots
subplot(3,1,1);
plot(Z);
    title('Low-frequency input');
subplot(3,1,2);
stem(x,'r');
    title ('Raw interpolation (padding with zeros)');
subplot(3,1,3);
plot([z w]); 
    legend('final','intermediate','Location','best'); 
    title('Low-pass interpolation')
