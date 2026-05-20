function Z = desvec(Zv,M)
% PURPOSE: Creates a matrix unstacking a vector
% ------------------------------------------------------------
% SYNTAX:  Z = desvec(Zv,M);
% ------------------------------------------------------------
% OUTPUT:  Z = matrix with M cols.
% ------------------------------------------------------------
% INPUT: Zv    : a vector nx1 
%	     M	   : number of cols.
% ------------------------------------------------------------

% written by:
%  Enrique M. Quilis

% Version 1.1 [December 2018]

[n,m] = size(Zv);

if (m ~= 1)
   error (' *** NUMBER OF COLUMNS GREATER THAN ONE *** ')
end

N = round(n/M);

Z(:,1) = Zv(1:N,1);
j=2;
while (j <= M)
   Z = [ Z  Zv( (j-1)*N+1 : j*N , 1)];
   j=j+1;
end
