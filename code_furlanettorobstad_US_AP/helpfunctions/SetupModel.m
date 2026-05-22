function setup = SetupModel(num)

% user choices
setup.nlag         = 12;                   % number of lags
setup.nex          = 1;                    % constant term
setup.horizon      = 60;                   % maximum horizon for IRFs
setup.horizons     = [0,1,2,3,4,inf];     % horizons upon which sign and zero restrictions can be imposed
setup.M0           = 1e4;                  % number of draws for Gibbs sampler
setup.iter_show    = 1e2;                  % display option to assess the number of draws obtained
setup.prior_type   = 'asymmetric';         % 'flat' (i.e., uhlig 2005) or 'minnesota'
setup.label_R      = 'baseline';           % label for selected identification scheme
setup.reoptimizeMH = 1;                    % set to 1 to reoptimize hyper parameters of normal-inverse wishart prior
setup.save_every   = 10;                   % thinning, we save every other "save_every" draw
setup.nsave        = setup.M0/setup.save_every;       % we store nsave elements
setup.elliptical_Q = 1;                    % set to  0 to use Chan, Matthes, Yu instead of ESS for Q.
setup.seed         = 0;                    % reproduce seed from housekeeping. Useful for exact replication across different OS
setup.nshocks      = 4;                   % number of shocks to identify, up to 10. To add more than 10 see identification section
setup.load_existing_init = 1;             % load existing initial values
% endogenous settings
setup.nvar   = size(num,2);                       % number of endogenous variables
setup.m      = setup.nvar*setup.nlag + setup.nex; % number of exogenous variables
setup.NS     = 1 + numel(setup.horizons);       % number of objects in F(L_0,L_+) to which we impose sign restrictions: F(L_0,L_+)=L_0;
setup.e      = eye(setup.nvar);                 % identity matrix
setup.h      = @(x)chol(x);                     % upper triangular Cholesky factorization of \Sigma



%==========================================================================
%% identification: 
% L_0 is the matrix of contemporaneous impulse responses
% L_0(i,j) is the response of variable i to shock j
% L_0(i,j) S_{i,j} L_0 e_j 
%==========================================================================
setup.Ss = cell([setup.nvar,1]);
for ii=1:setup.nvar
    setup.Ss{ii}=zeros(0,setup.nvar*setup.NS);
end


switch setup.label_R
    case 'baseline'
        baseline_identification;
end




