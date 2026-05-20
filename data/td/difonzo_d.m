% PURPOSE: Demo of difonzo()
%          Temporal disaggregation with indicators.
%          Multivariate model with transversal constraint
% 			  di Fonzo method
%---------------------------------------------------
% USAGE: dinfonzo_d
%--------------------------------------------------

close all; clear all; clc;

% Loading data
load('Spain_employment_regional_4.mat')

% ---------------------------------------------
% Inputs

% Type of aggregation
ta = 2;
% Frequency conversion 
sc = 4;    
% Model for the innovations: white noise (0), random walk (1)
type = 1;
% Number of high frequency indicators linked to each low frequency
% aggregate
f = ones(1,size(x,2));
% Multivariate temporal disaggregation
res = difonzo(Y,x,z,ta,sc,type,f);
% Printed output
file_out = 'difonzo.out';   
tdprint(res,file_out);
edit difonzo.out;
% Graphs
tdplot(res);

figure
vnames = {'NORTH NORTHWEST', 'SOUTH SOUTHEAST', 'CENTER', 'ISLANDS'};
for j=1:size(x,2)
    subplot(2,2,j)
    plot([res.y(:,j)],'-b');
    hold on
    plot([res.y(:,j)-res.d_y(:,j) res.y(:,j)+res.d_y(:,j)],':r');
    title(vnames(j));
end

