function ztab = tab_s(z,s)
% PURPOSE: Write a time series in matrix form: years x seasons
% ------------------------------------------------------------
% SYNTAX: ztab = tab_s(z,s); 
% ------------------------------------------------------------
% OUTPUT: ztab: Nxs -> matrix form: years x seasons
% ------------------------------------------------------------
% INPUT:  z: nx1 -> time series
%         s: 1x1 -> seasonal frequency
% ------------------------------------------------------------
% SEE ALSO: acc, de_acc
% ------------------------------------------------------------
% NOTE: This function assumes that the first year is complete. The last
% year can be incomplete.

% written by:
%  Enrique M. Quilis

% Version 1.0 [December 2018]

% Data dimension
n = length(z);
N = ceil(n/s);

% Implicit data dimension (if the last year is complete)
m = N*s;

% Completing matrix if the last year is not complete
if (m > n)
    z(end+1:m) = NaN;
end

% Forming matrix: rows = years, cols. = months/quarters
ztab = (desvec(z,N))';