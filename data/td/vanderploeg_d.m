% PURPOSE: Demo of vanderploeg()
%          Estimation subject to transversal constraints
%          by means of a quadratic optimization criterion
%-----------------------------------------------------------------
% USAGE: vanderploeg_d
%-----------------------------------------------------------------

close all; clear all; clc;

%-----------------------------------------------------------------
% Unbalanced cross-section vector
y = [   220.00
    130.00
    200.00
    100.00
    450.00
    70.00
    120.00
    221.00 ];

[k, n]=size(y);

%-----------------------------------------------------------------
% Linear constraints
A =[     1.00             0
    1.00             0
    1.00          1.00
    1.00             0
    -1.00             0
    -1.00             0
    -1.00             0
    -1.00         -1.00  ];

[k, m] = size(A);

%-----------------------------------------------------------------
% VCV matrix of estimates

% Vector of variances
% Note: Fixed estimation: s(5)=0 --> z(5)=y(5)
s = [10 5 25 55 0 15 10 12];
Aux1 = (diag(sqrt(s)));

% Correlation matrix: C
C = zeros(k);
C(1,3) = 0.5;
Aux2 = tril(C');
C = C + Aux2 + diag(ones(1,k));

% VCV matrix: S
S = Aux1 * C * Aux1;

% van der Ploeg balancing
res = vanderploeg(y,S,A);

% Check
format bank
disp (''); disp ('*** INITIAL AND FINAL DISCREPANCIES ***'); disp('');
[ A' * y  A' * res.z]

% Revision (as %)
p = 100 * ((res.z - y) ./ y);

% Final results:
disp ('');
disp ('*** INITIAL ESTIMATE, FINAL ESTIMATE, REVISION AS %, INITIAL VARIANCES, FINAL VARIANCES ***');
disp ('');
[y res.z p diag(S) diag(res.Sz)]
format short

% Graphs
sv = (diag(res.Sz));
s  = diag(S);
for j=1:k
    if (s(j) == 0)
        x = linspace(min(y(j),res.z(j))*0.9,max(y(j),res.z(j))*1.1,1000);
        subplot(4,2,j)
        plot(x,0)
        hold on
        plot(y(j),0,'or','LineWidth',6)
        axis([min(y(j),res.z(j))*0.9 max(y(j),res.z(j))*1.1 0 0.2])
        legend('Absolute tight prior','Location','best')
    else
        x = linspace(min(y(j),res.z(j))*0.8,max(y(j),res.z(j))*1.2,1000);
        y0 = (1/ (sqrt(2*pi*s(j)))) * exp(-((x-y(j)).^2 ./ s(j)));
        y1 = (1/ (sqrt(2*pi*sv(j)))) * exp(-((x-res.z(j)).^2 ./ sv(j)));
        subplot(4,2,j)
        plot(x, y0, x, y1,'-r')
        if j==2
        legend('Prior','Posterior','Location','best')
        end
    end
end