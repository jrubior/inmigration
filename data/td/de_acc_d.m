% PURPOSE: demo of de_acc()
%          Extracting original time series from accumulated time series
%---------------------------------------------------
% USAGE: de_acc_d
%---------------------------------------------------

close all; clear all; clc;

% Loading data
load bournay_laroque;
za = x;

% Seasonal frequency
s = 4;

% Calling function
z = de_acc(za,s);

% Plot
plot([za z]);
legend('Accumulated','De-accumulated','Location','best');
