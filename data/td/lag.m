function X = lag(Z,p)
% PURPOSE: Create a matrix with the lags of Z up to p
% ------------------------------------------------------------
% SYNTAX: X = lag(Z,p)
% ------------------------------------------------------------
% OUTPUT: X = matrix of Z(1)..Z(p), including a vector of 1s,
%				 in the common range of obs.: p+1 .. n
%			    X : n-p x k*p+1
% ------------------------------------------------------------
% INPUT: Z = an nxk matrix of original series, columnwise
%			p = the maximum lag, p>=1

% ------------------------------------------------------------
% written by:
%  Enrique M. Quilis

% Cheks
if (p < 1)
   error (' *** The number of lags should be greater than one *** ' );
end

% Data dimension
[n,k] = size(Z);

if (n < p)
   error (' *** The # of obs. should be greater than the # of lags *** ' );
end

% Initial matrix
XX = [Z
   zeros(p,k)];

% Recurrence
j = 1;
while (j <= p)
   Aux = [zeros(j,k)
            Z
   		  zeros(p-j,k)];
   XX = [ XX Aux ];
   j = j + 1;
end

% Final matrix
Xf = [zeros(p,k)
     Z];
  
% Adding the last lag

XX = [XX Xf];

% Selecting the common range

X = XX(p+1:n,k+1:k*(p+1));

% Preparing X to include a constant

[nX,mX] = size(X);

X = [ ones(nX,1) X];

% end of function lag
