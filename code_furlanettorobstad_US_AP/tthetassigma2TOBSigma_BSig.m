function [Bdraw,Sigmadraw,A,Bstrdraw] = tthetassigma2TOBSigma_BSig(ttheta_old,ssigma_old,info,ii,Bold,Sigold,Aold,Bstrold)

% Assuming you have ttheta_old (cell array) and ssigma_old (vector) as inputs
% and nvar is defined

nvar = info.nvar;
% Initialize matrices
Bdraw = Bstrold;
A = Aold;

% Extract B coefficients and alpha parameters from ttheta_old
count_alp = 0;
i = ii;
% for i = 1:nvar
    if i == 1
        % First equation: only B coefficients (no alpha terms)
        Bdraw(:,i) = ttheta_old{i};
    else
        % Other equations: B coefficients + alpha terms
        n_alpha = i - 1;  % Number of alpha terms for equation i
        n_beta = length(ttheta_old{i}) - n_alpha;  % Number of B coefficients
        
        % Extract B coefficients (first part)
        Bdraw(:,i) = ttheta_old{i}(1:n_beta);
        
        % Extract alpha coefficients (last part)
        aalpha_i = ttheta_old{i}(n_beta+1:end);
        
        % Place alpha coefficients in matrix A
        A(i, 1:i-1) = aalpha_i';
    end
% end

% Reconstruct Sigmadraw
% Create diagonal matrix D from ssigma_old
D = diag(ssigma_old);

% Get A inverse
Ainv = A \ eye(nvar);

% Create lower triangular matrix C
C = Ainv * sqrt(D);

% Reconstruct Sigmadraw using C*C'
Sigmadraw = C * C';

Bstrdraw = Bdraw; %structural
%Bdraw = Bdraw/(A'); %***we can reuse Ainv, reduced form
Bdraw = Bdraw*Ainv';