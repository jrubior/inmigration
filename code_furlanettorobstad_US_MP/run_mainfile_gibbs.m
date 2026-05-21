%==========================================================================
%% housekeeping
%==========================================================================
clear variables;close all;userpath('clear');restoredefaultpath;clc;
tic;
rng('default');                                % reinitialize the random number generator to its startup configuration
seed =0; rng(seed,'twister');                  % set seed
currdir=pwd;
addpath([currdir,'/data']);                    % set path to data functions
addpath(genpath([currdir,'/helpfunctions/'])); % set path to helper functions
addpath('ess_within_smc');                     % set path to helper functions to find initial values for the Gibbs Sampler
cd(currdir)


%==========================================================================
%% load the data
%==========================================================================

% Get the folder where this m-file lives
scriptDir = fileparts(mfilename('fullpath'));

% Build path to the data file (go up one level, then into 'data')
dataFile = fullfile(scriptDir, '..', 'data', 'data_FURLANETTO_ROBSTAD_USA.csv');

% Read the CSV
T    = readtable(dataFile);
data = table2array(T);

var_id = 1:5;             % index for variables included in the SVAR
idx_ns = 1:5;             % index for variables in levels   
num    = data(:,var_id);  % data that enters the estimation

%==========================================================================
%% model setup
%==========================================================================
setup = SetupModel(num);

%==========================================================================
%% mappings
%==========================================================================
fo                 = @(x)f_h(x,setup);              % (vec(A_0),vec(A_+))        -> (vec(B),vec(Sigma),vec(Q))
fo_inv             = @(x)f_h_inv(x,setup);          % (vec(B),vec(Sigma),vec(Q)) -> (vec(A_0),vec(A_+))
fo_str2irfs        = @(x)StructuralToIRF(x,setup);  % (vec(A_0),vec(A_+))        -> (vec(L_0),vec(L_+))
fo_str2irfs_inv    = @(x)IRFToStructural(x,setup);   % (vec(L_0),vec(L_+))        -> ((vec(A_0),vec(A_+))
gs_qr              = @(x)qr_unique(x);              % (X)                        -> Q, where X(i,j) is standard normal and Q is an orthogonal matrix. Stewart (1980) QR factorization

%==========================================================================
%% model representation
%==========================================================================
% Rubio, Waggoner, and Zha (RES 2010)'s notation
% yt(t) A0 = xt(t) Aplus + constant + et(t) for t=1...,T;
% yt(t)    = xt(t) B     + ut(t)            for t=1...,T;
% x(t)     = [constant,yt(t-1), ... ,yt(t-nlag)];
% matrix notation yt = xt*B + ut;
% xt=[yt_{-1} ones(T,1)];
% reduced-form representation: Y=X*B + U
% Y = [yt(1);...;yt(T)]
% X = [xt(1);...;xt(T)]

Y  = num(setup.nlag+1:end,var_id);
[T,nvar] = size(Y);
X = zeros(T,setup.nvar*setup.nlag); 
for ii=1:setup.nlag
    X(:,(ii-1)*setup.nvar+1:ii*setup.nvar) = num(setup.nlag-ii+1:end-ii,:);
end
X = [ones(T,1) X];


% when setting hyperparameters sometimes dummy observations are stack to observations, hence we use Ynd, Xnd to refer to the data before dummy observations are introduced are added.
Ynd = Y;
Xnd = X;

%==========================================================================
%% Prior 
%==========================================================================
% Prior for reduced-form parameters (B,Sigma) is Normal-Inverse-Wishart, prior for Q is uniform under Haar Measure
switch setup.prior_type

    case 'flat' % Uhlig (2005): bound to fail in large SVARs

        %% prior
        nnuBar             = 0;
        PphiBar            = zeros(setup.nvar);   
        PpsiBar            = zeros(setup.m,setup.nvar);  
        OomegaBarInverse   = zeros(setup.m);

    case 'minnesota'
        
        [nnuBar,PphiBar,PpsiBar,OomegaBarInverse,Y,X] = SetupMinnesota(Ynd,Xnd,setup);

end



%==========================================================================
%% Posterior
%==========================================================================
% posterior for reduced-form parameters
setup.nnuTilde            = T +nnuBar;
setup.OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(setup.m);
setup.OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
setup.PpsiTilde           = (X'*X  + OomegaBarInverse)\(X'*Y + OomegaBarInverse*PpsiBar);
setup.PphiTilde           = Y'*Y + PphiBar + PpsiBar'*OomegaBarInverse*PpsiBar - setup.PpsiTilde'*setup.OomegaTildeInverse*setup.PpsiTilde;
setup.PphiTilde           = (setup.PphiTilde'+setup.PphiTilde)*0.5;




%==========================================================================
%% Set objects to store draws and useful definitions
%==========================================================================
% definitios used to store orthogonal-reduced-form draws
Bdraws         = cell([setup.nsave,1]); % reduced-form lag parameters
Sigmadraws     = cell([setup.nsave,1]); % reduced-form covariance matrices
Qdraws         = cell([setup.nsave,1]); % orthogonal matrices


% definitions to facilitate the draws from Sigma
Rmean = nan(nvar,setup.nnuTilde);
for i=1:setup.nnuTilde
    Rmean(:,i) = zeros(setup.nvar,1);
end
Rvariance = zeros(setup.nvar*setup.nnuTilde,setup.nvar*setup.nnuTilde);
for i=1:setup.nnuTilde
    Rvariance((i-1)*setup.nvar+1:i*setup.nvar,(i-1)*nvar+1:i*nvar) = setup.PphiTilde\eye(setup.nvar);
end
chol_lower_Rvariance = chol(Rvariance,'lower');
chol_lower_Rvariance_common = chol(setup.PphiTilde\eye(setup.nvar), 'lower');

% definitions to facilitate the draws from B|Sigma
hh              = setup.h;
cholOomegaTilde = hh(setup.OomegaTilde)'; % this matrix is used to draw B|Sigma below
chol_lower_OomegaTilde     = chol(setup.OomegaTilde, 'lower');
chol_lower_OomegaTilde_inv = eye(size(setup.OomegaTilde))/chol_lower_OomegaTilde;


% definitions related to IRFs and stability of the coefficients; based on page 12 of Rubio, Waggoner, and Zha (RES 2010)
setup.J      = [setup.e;repmat(zeros(setup.nvar),setup.nlag-1,1)];
setup.A      = cell(setup.nlag,1);
setup.extraF = repmat(zeros(setup.nvar),1,setup.nlag-1);
setup.bigF   = zeros(setup.nlag*setup.nvar,setup.nlag*setup.nvar);
for l=1:setup.nlag-1
    setup.bigF((l-1)*setup.nvar+1:l*setup.nvar,setup.nvar+1:setup.nlag*setup.nvar)=[repmat(zeros(setup.nvar),1,l-1) setup.e repmat(zeros(setup.nvar),1,setup.nlag-(l+1))];
end

% define restrictions to be used when truncating likelihood
switch setup.label_R
    case 'baseline'
        function_restrictions = @restriction_baseline;
end

%==========================================================================
%% Gibbs Initialization
%==========================================================================
if setup.load_existing_init ==1
    load(['init',setup.label_R,'.mat'],'z_old_B','z_old_R','z_old_X','Bdraw','Sigmadraw','Qdraw')
else
disp('initilization of Q ...');
init_time =tic;
function_restrictions_Svec = @SRvec;
rst = ess_within_smc_init(num, function_restrictions_Svec,setup);
z_old_B   = vec(squeeze(rst.mat_B(:,:,1)));
z_old_R   = rst.mat_R(:,1);
z_old_X   = rst.mat_X(:,1);
Bdraw     = rst.mat_B(:,:,1);
Sigmadraw = rst.mat_S(:,:,1);
Qdraw = rst.mat_Q(:,:,1);
    save(['init',setup.label_R,'.mat'],'z_old_B','z_old_R','z_old_X','Bdraw','Sigmadraw','Qdraw')
disp('initilization of Q ... done ...');
time_to_initialize = toc(init_time);
disp(time_to_initialize)
end



%% initialize counters to track the state of the computations
counter = 1; % counter to display
record  = 1; % counter of draws


% check initial values satisfy restrictions
function_restrictions_i = @(ag0, ag1, ag2, ag3) function_restrictions(ag0, ag1, ag2, ag3);
S_prop = function_restrictions_i(Bdraw, Sigmadraw, Qdraw,setup);
if S_prop==0
disp('Initial values are not valid')
return;
end


tstart = tic;

ind_save=0;

while record<=setup.M0





    %% Using Elliptical Slice Sampling, draw Q^i
    slice.scale_z    = 1;
    slice.mean       = zeros(setup.nvar*setup.nvar,1);
    slice.chol_cov_z = chol(eye(setup.nvar*setup.nvar),'lower');
    slice.chol_cov_z = 1;
    slice.fcn_lik    = @(z_prop) loglike_Q(z_prop,function_restrictions_i,gs_qr,Bdraw,Sigmadraw,setup);
    slice.nobs       = setup.nvar*setup.nvar;
    % --- actual slice sampling
    lik_old = slice.fcn_lik(z_old_X);
    [z_old_X, ~, n_try_Q] = slice_sampling_v02(slice, z_old_X, lik_old);
    [Qdraw,Rdraw] =  gs_qr(reshape(z_old_X,setup.nvar,setup.nvar));



    %% Using Elliptical Slice Sampling, draw Sigma
    slice.scale_z    = 1;
    slice.mean       = vec(Rmean);
    slice.chol_cov_z_common = chol_lower_Rvariance_common;slice.nobs1 = setup.nvar;
    slice.nobs2      = setup.nnuTilde;
    slice.fcn_lik    = @(z_prop) loglike_Sigma_mc_v02(z_prop,function_restrictions_i,Bdraw,Qdraw,chol_lower_OomegaTilde,chol_lower_OomegaTilde_inv,setup);
    slice.nobs       = setup.nvar*setup.nnuTilde;
    % --- actual slice sampling
    lik_old = slice.fcn_lik(z_old_R);
    [z_old_R, ~, n_try_Sigma] = slice_sampling_v02_Sigma(slice, z_old_R, lik_old);
    R_old = reshape(z_old_R,setup.nvar,setup.nnuTilde);
    Sigmadraw=inv(R_old*R_old');


    %% Using Elliptical Slice Sampling, draw B
    cholLSigmadraw = setup.h(Sigmadraw)';
    slice.scale_z    = 1;
    slice.mean       = vec(setup.PpsiTilde);
    slice.chol_cov_z1 = cholLSigmadraw;
    slice.chol_cov_z2 = cholOomegaTilde;
    slice.nobs1 = setup.nvar;
    slice.nobs2 = setup.m;
    slice.fcn_lik    = @(z_prop) loglike_B(z_prop,function_restrictions_i,Sigmadraw,Qdraw,setup);
    slice.nobs       = setup.nvar*setup.m;
    % --- actual slice sampling
    lik_old = slice.fcn_lik(z_old_B);
    [z_old_B, ~, n_try_B] = slice_sampling_v02_kron(slice, z_old_B, lik_old);
    Bdraw = reshape(z_old_B,setup.m,setup.nvar);


    if rem(record, setup.save_every) == 0
        % --- store
        % store orthogonal reduced-form draws
        ind_save = ind_save + 1;
        Bdraws{ind_save,1}     = Bdraw;
        Sigmadraws{ind_save,1} = Sigmadraw;
        Qdraws{ind_save,1} = Qdraw;
    end

    record=record+1;

    if counter==setup.iter_show

        display(['Number of draws = ',num2str(record)])
        display(['Remaining draws = ',num2str(setup.M0-(record))])
        counter =0;

    end
    counter = counter + 1;
end

telapsed = toc(tstart);

disp(telapsed)

%% compute and store impulse responses (which are our objects of interest)
L    = zeros(setup.horizon+1,setup.nvar,setup.nshocks,ind_save);
for s=1:ind_save
    Bdraw =     Bdraws{s,1} ;
    Sigmadraw = Sigmadraws{s,1} ;
    Qdraw=Qdraws{s,1};
    BSigmaQ = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];
    structpara = f_h_inv(BSigmaQ,setup);
    LIRF = IRF_horizons(structpara, setup.nvar, setup.nlag, setup.m, setup.nex, 0:setup.horizon);
    for h=0:setup.horizon
        L(h+1,:,1:setup.nshocks,s) =  LIRF(1+h*nvar:(h+1)*setup.nvar,1:setup.nshocks);
    end
end



cd results

if ~exist('matfiles', 'dir')
    mkdir('matfiles')
end

savefile= ['matfiles/rgibbs_',setup.label_R,'prior_',setup.prior_type,'_ndraws',num2str(setup.M0),'number_of_shocks',num2str(setup.nshocks),'save_every',num2str(setup.save_every),'.mat'];

save(savefile,'L','setup','telapsed','-v7.3');
cd ..
