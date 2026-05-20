% PURPOSE: demo of tab_s()
%          Matrix form of a time series. Each column is a month/quarter.
%---------------------------------------------------
% USAGE: tabs_s_d
%---------------------------------------------------

close all; clear all; clc;

% Loading data
load bournay_laroque;
z = x;
     
% Seasonal frequency
s = 4;

% Calling function
z_tab = tab_s(z,s);

% Display
disp('Rows are years, cols. are months');
disp('');
z_tab
