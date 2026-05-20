function dz = de_acc(z,s)
% PURPOSE: De-accumulate a YTD accumulated time series
% ------------------------------------------------------------
% SYNTAX: dz = de_acc(z,s);
% ------------------------------------------------------------
% OUTPUT: dz: nx1 -> de-accumulated time series
% ------------------------------------------------------------
% INPUT:  z: nx1 -> accumulated time series
%         s: 1x1 -> seasonal frequency
% ------------------------------------------------------------
% LIBRARY: tab_s
% ------------------------------------------------------------
% SEE ALSO: acc, accYTD
% ------------------------------------------------------------

% written by:
%  Enrique M. Quilis

% Version 1.0 [December 2018]

% Forming matrix: rows = years, cols. = months/quarters
ztab = tab_s(z,s);

% Differencing the matrix columnwise generated de de-accumulated time
% series, still in matrix form
dZ = ztab(:,1);
for j=2:s
    dZ(:,j) = ztab(:,j) - ztab(:,j-1);
end

% Transforming the matrix into a vector
dz = vec(dZ');

% Deleting NaN observations that may appear at the end of the sample 
I = isnan(dz);
dz(I==1) = [];

