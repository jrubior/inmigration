function [g1 g2] = ar_order(z,pmax,opG)
% PURPOSE: Determination of AR order through AIC and BIC
% ------------------------------------------------------------
% SYNTAX:  [g1 g2] = ar_order(z,pmax,opG);
% ------------------------------------------------------------
% OUTPUT: g1: pmax x 1 -> AIC
%         g2: pmax x 1 -> BIC
%         Two plots: AIC vs order, BIC vs order
% ------------------------------------------------------------
% INPUT: z        = an nx1 matrix original series (columnwise)
%        pmax     = maximum order of AR model
%        opG: 1x1 -> Optional: graphics: yes (1), no (0). Default = 0.
% ------------------------------------------------------------
% LIBRARY: arx
% ------------------------------------------------------------
% SEE ALSO: order_var

% ------------------------------------------------------------
% written by:
%  Enrique M. Quilis
%
% December 2011
% Version 1.2

% ------------------------------------------------------------
% Setting default: opG=0 (no graphical output)
if (nargin == 2)
    opG = 0;
end

% Initial allocation
g1 = zeros(pmax,1);
g2 = zeros(pmax,1);

% ------------------------------------------------------------
% Checks

[n,k]   = size(z);
if (k > 1)
    error('*** ONLY UNIVARIATE MODELING, PLEASE ***');
end

% ------------------------------------------------------------
% Loading the structure

rex.z = z;
rex.X = [];
rex.opC = 1;

% ------------------------------------------------------------
% Basic loop

for p=1:pmax
   rex.p = p;
   res = arx(rex);
   g1(p) = res.perform(2);
   g2(p) = res.perform(3);
end

% ------------------------------------------------------------
% Optional graphics

switch opG
case 1
    subplot(1,2,1); plot(g1,'-bs'); title('AIC'); xlabel('AR order'); grid on;
    subplot(1,2,2); plot(g2,'-rs'); title('BIC'); xlabel('AR order'); grid on;
end

