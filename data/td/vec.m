function Y_big = vec(Y)
% PURPOSE: Create a matrix stacking the columns of Y
% ------------------------------------------------------------
% SYNTAX: Y_big = vec(Y);
% ------------------------------------------------------------
% OUTPUT: Y_big = matrix of columns of Y
% ------------------------------------------------------------
% INPUT: Y = an nxM matrix of original series, columnwise
% ------------------------------------------------------------

% NOTE: An alternative and more efficient version is available from James
% LeSage "Econometrics Toolbox".

% written by:
%  Enrique M. Quilis

% Version 1.1 [December 2018]

% Data dimension
M = size(Y,2);

Y_big = Y(:,1);   % Starting the recursion
j = 2;
while (j <=  M)
   Y_big = [Y_big
          Y(:,j) ];
   j = j+1;
end