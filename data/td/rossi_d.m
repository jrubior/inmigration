% PURPOSE: Demo of rossi()
%          Temporal disaggregation with indicators.
%          Multivariate model with transversal constraint
% 			  Rossi method
%---------------------------------------------------
% USAGE: rossi_d
%---------------------------------------------------

close all; clear all; clc;

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
% Type of univariate disaggregation procedure
opMethod = 2;
% Type of univariate disaggregation procedure: estimation method
type = 1;
% Multivariate temporal disaggregation
res = rossi(Y,x,z,ta,sc,opMethod,type);
% Printed output
file_out = 'rossi.out';   
tdprint(res,file_out);
edit rossi.out;
% Graphs
tdplot(res);
figure
vnames = {'NORTH NORTHWEST', 'SOUTH SOUTHEAST', 'CENTER', 'ISLANDS'};
for j=1:size(x,2)
    subplot(2,2,j)
    plot([x(:,j) res.y_prelim(:,j) res.y(:,j)]);
    title(vnames(j));
    legend('HF Tracker','Prelim. Estimate','Estimate','Location','best');
end



