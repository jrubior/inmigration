% PURPOSE: Demo of bal()
%          Proportional transversal balancing
%---------------------------------------------------
% USAGE: bal_d
%---------------------------------------------------

clc; clear all; close all;

% Unbalanced time series (read columnwise)
y = [120    20      55
    130    22      66
    190    33      99
    170    35      110 ];

% Transversal constraints
z=[ 210
    208
    300
    310 ];

% Balancing
yb = bal(y,z);

% Graphs
t = (1:size(y,1));
subplot(2,1,1)
plot(t, y, '-b', t, yb, '-r')
    xlabel('Time')
    title('Unbalanced (blue) and balanced (red) time series')
subplot(2,1,2)
plot(t, sum(y')', '-b', t, z, '-r')
    xlabel('Time')
    legend('Unbalanced total', 'Constraint (red)','Location','best')
    