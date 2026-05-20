% PURPOSE: demo of auxiliary functions
%-----------------------------------------------------------------------
% USAGE: aggreg_d
%-----------------------------------------------------------------------

close all; clear all; clc;

% Loading data
load bournay_laroque;

%====================================================================================
%  TEMPORAL AGGREGATION: MONTHLY --> ANNUAL
%====================================================================================

% Type of aggregation
op1 = 1;
% Frequency conversion
sc = 4;
% Calling function
X = temporal_agg(x,op1,sc);


% Graphs
subplot(2,1,1);
t = 1:length(x);
plot(t,x,'or-');
    title('Original quarterly data');
subplot(2,1,2);
t = 1:length(X);
plot(t,X,'ob-');
    title('Annual aggregation');

%====================================================================================
%  NAIVE TEMPORAL DISAGGREGATION: ANNUAL --> MONTHLY
%====================================================================================

y = copylow(Y,op1,sc);
figure;
% Graphs
subplot(2,1,1);
t = 1:length(Y);
plot(t,Y,'ob-');
    title('Original annual data');
subplot(2,1,2);
t = 1:length(y);
plot(t,y,'or-');
    title('Quarterly copy');

% Inserting annual data in a string of zeros
y = copylow(Y,3,sc);
figure;
% Graphs
subplot(2,1,1);
t = 1:length(Y);
plot(t,Y,'ob-');
    title('Original annual data');
subplot(2,1,2);
t = 1:length(y);
plot(t,y,'or-');
    title('Inserting annual data in a string of zeros');

%====================================================================================
%    SYSTEMATIC SAMPLING OF HIGH FREQUENCY DATA
%====================================================================================

xs = ssampler(x,op1,sc);
figure;
% Graphs
subplot(2,1,1);
t = 1:length(x);
plot(t,x,'or-');
    title('Systematic sampling');
subplot(2,1,2);
t = 1:length(xs);
plot(t,xs,'ob-');
    title('');

%====================================================================================
%    YTD ACCUMULATION OF HIGH FREQUENCY DATA
%====================================================================================

xa = temporal_acc(x,op1,sc);
figure;
% Graphs
subplot(2,1,1);
t = 1:length(x);
plot(t,x,'or-');
    title('YTD Accumulation');
subplot(2,1,2);
t = 1:length(xa);
plot(t,xa,'ob-');
    title('');
