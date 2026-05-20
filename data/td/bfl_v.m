function y = bfl_v(Y,ta,d,sc)
% PURPOSE: Temporal disaggregation of a vector time series using the BFL method
% See bfl() for details. Note: the inputs are common across series.

% written by:
%  Enrique M. Quilis

% Version 1. [December 2018]

% Data dimension
[N, k] = size(Y);

% Cross-section loop
y = nan(N*sc, k);
for j=1:k
    res = bfl(Y(:,j),ta,d,sc);
    y(:,j) = res.y;
end
