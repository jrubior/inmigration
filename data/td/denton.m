function res = denton(Y,x,ta,d,sc,op1)
% PURPOSE: Temporal disaggregation using the Denton method.
% -----------------------------------------------------------------------
% SYNTAX: res = denton(Y,x,ta,d,sc,op1,op2);
% -----------------------------------------------------------------------
% OUTPUT: res: a structure
%         res.meth  = 'Denton';
%         res.N     = Number of low frequency data
%         res.ta    = Type of disaggregation
%         res.sc    = Frequency conversion
%         res.d     = Degree of differencing
%         res.y     = High frequency estimate
%         res.x     = High frequency indicator
%         res.U     = Low frequency residuals
%         res.u     = High frequency residuals
%         res.et    = Elapsed time
% -----------------------------------------------------------------------
% INPUT: Y: Nx1 ---> vector of low frequency data
%        x: nx1 ---> vector of low frequency data
%        ta: type of disaggregation
%            ta=1 ---> sum (flow)
%            ta=2 ---> average (index)
%            ta=3 ---> last element (stock) ---> interpolation
%            ta=4 ---> first element (stock) ---> interpolation
%        d: objective function to be minimized: volatility of ...
%            d=1 ---> first differences
%            d=2 ---> second differences
%        sc: number of high frequency data points for each low frequency data point
%            Some examples:
%            sc= 4 ---> annual to quarterly
%            sc=12 ---> annual to monthly
%            sc= 3 ---> quarterly to monthly
%        op1: additive variant (1) or proportional variant(2)
% -----------------------------------------------------------------------
% LIBRARY: aggreg, dif
% -----------------------------------------------------------------------
% SEE ALSO: tdprint, tdplot
% -----------------------------------------------------------------------
% REFERENCE: Denton, F.T. (1971) "Adjustment of monthly or quarterly
% series to annual totals: an approach based on quadratic minimization",
% Journal of the American Statistical Society, vol. 66, n. 333, p. 99-102.

% ------------------------------------------------------------
% written by:
%  Enrique M. Quilis
% Version 3.0 [December 2018]

t0 = clock;

% -----------------------------------------------------------------------
% Checking d (minimization of the volatility of first differences or second
% differences).
if ((d ~= 1) & (d ~= 2))
    error ('*** d must be 1 or 2 ***');
end

% -----------------------------------------------------------------------
% Computing input dimensions.
[N,M] = size(Y);
[n,m] = size(x);

% Ckecking consistency among dimensions.
if ((m > 1) | (M > 1))
    error ('*** x HAS INCORRECT DIMENSION: col(x)>1 ***');
else
    clear m M;
end

% Defining combination of variant (additive or proportional) and degree of
% differencing (d=1 or d=2).
FLAX = [1 2; 3 4];
flax = FLAX(op1,d);

% -----------------------------------------------------------------------
% Generation of aggregation matrix C: Nxn.
C = aggreg(ta,N,sc);

% -----------------------------------------------------------
% Expanding the aggregation matrix to perform
% extrapolation if needed.
if (n > sc * N)
    pred = n - sc*N;           % Number of required extrapolations
    C = [C zeros(N,pred)];
else
    pred = 0;
end

% -----------------------------------------------------------------------
% Temporal aggregation matrix of the indicator: Nx1.
X = C*x;

% -----------------------------------------------------------------------
% Computing low frequency residuals: Nx1.
U = Y - X;

% -----------------------------------------------------------------------
% Computing difference matrix: nxn. (Note: initial conditions = 0).
D = dif(1,n);

switch flax
    case 1 %Additive first differences: d=1.
        Q = inv(D'*D);
    case 2 %Additive second differences: d=2.
        D2 = dif(2,n);
        Q = inv(D2'*D2);
    case 3 %Proportional first differences: d=1.
        Q = diag(x) * inv(D'*D) * diag(x);
    case 4 %Proportional second differences: d=2.
        D2 = dif(2,n);
        Q = diag(x) * inv(D2'*D2) * diag(x);
end %of switch flax.

% High frequency residuals.
u = Q*C' * inv(C*Q*C') * U;

% u = Qi*C'*inv(C*Qi*C')*U;
% High frequency estimator.
y = x + u;

% -----------------------------------------------------------------------
% Loading the structure
% -----------------------------------------------------------------------
% Basic parameters
res.meth = 'Denton';
switch op1
    case 1
        res.meth1 = 'Additive variant';
    case 2
        res.meth1 = 'Proportional variant';
end
res.N 	= N;
res.ta	= ta;
res.sc 	= sc;
res.pred = pred;
res.d 	= d;
% -----------------------------------------------------------------------
% Series
res.y = y;
res.x = x;
res.U = U;
res.u = u;
% -----------------------------------------------------------------------
% Elapsed time
res.et  = etime(clock,t0);
