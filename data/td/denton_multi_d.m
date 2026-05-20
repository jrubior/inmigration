% PURPOSE: Demo of denton()
%          Temporal disaggregation with indicators.
%          Multivariate model with transversal constraint
% 		   Denton method, addititive or proportional variants.
%---------------------------------------------------
% USAGE: denton_d
%---------------------------------------------------

clc; close all; clear all;

% Loading data
load('Spain_employment_regional_4.mat')
x(end-2:end,:) = [];    %No free extrapolations
z(end-2:end) = [];      %No free extrapolations

% ---------------------------------------------
% Inputs

% Type of aggregation
ta = 2;
% Frequency conversion 
sc = 4;    
% Minimizing the volatility of d-differenced series.
d = 2;
% Additive (1) or proportional (2) variant [optional, default=1]
op1 = 1;
% Multivariate temporal disaggregation
res = denton_multi(Y,x,z,ta,sc,d,op1);
% Printed output
file_out = 'denton_multi.out';   
tdprint(res,file_out);
edit denton_multi.out;
% Graphs
tdplot(res);
figure
vnames = {'NORTH NORTHWEST', 'SOUTH SOUTHEAST', 'CENTER', 'ISLANDS'};
for j=1:size(x,2)
    subplot(2,2,j)
    plot([x(:,j) res.y(:,j) copylow(Y(:,j),1,sc)]);
    title(vnames(j));
    legend('HF Tracker','Estimate','LF Benchmark','Location','best');
end

