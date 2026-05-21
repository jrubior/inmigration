function Y = h_tilde(X, setup)
%
% X - n x n matrix
%
% Y - n x n matrix such that Y'*Y = 0.5*(X+X')
%

Y=setup.h(0.5*(X+X'));
