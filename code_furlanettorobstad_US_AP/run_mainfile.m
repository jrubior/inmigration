%==========================================================================
%% housekeeping
%==========================================================================
clear variables;close all;userpath('clear');restoredefaultpath;clc;
tic;
rng('default'); % reinitialize the random number generator to its startup configuration
seed =0;
rng(seed,'twister');         % set seed
currdir=pwd;
addpath([currdir,'/helpfunctions']); % set path to helper functions
addpath([currdir,'/helpfunctions/ChrisSimsOptimize']);
addpath([currdir,'/helpfunctions/subroutines']); % set path to helper functions
addpath([currdir,'/utility']); % implements Chan's asymmetric priors
cd(currdir)


%==========================================================================
%% load the data
%==========================================================================

% Get the folder where this m-file lives
scriptDir = fileparts(mfilename('fullpath'));

% Build path to the data file (go up one level, then into 'data')
dataFile = fullfile(scriptDir, '..', 'data', 'data_FURLANETTO_ROBSTAD_USA.csv');

% Read the CSV
T = readtable(dataFile);
data = table2array(T);

var_id = 1:5;             % index for variables included in the SVAR
idx_ns = 1:5;             % index for variables in levels   
num    = data(:,var_id);  % data that enters the estimation




%=========================================================================
%% model setup
%==========================================================================
nlag      = 12;                    % number of lags
nvar      = size(num,2);           % number of endogenous variables
nex       = 1;                     % 12 because of seasonal dummies
m         = nvar*nlag + nex;       % number of exogenous variables
horizon   = 60;                    % maximum horizon for IRFs
horizons  = 0;                     % horizons upon which sign and zero restrictions can be imposed
NS        = 1 + numel(horizons);   % number of objects in F(THETA) to which we impose sign and zero restrictions: F(THETA)=[A_0;L_{0};L_{1};L_{2};L_{3};L_{4};L_{inf}]
e         = eye(nvar);             % create identity matrix
M0        = 1e4;

iter_show  = 1e4;
label_R = 'cmy';

save_every   = 10;% save every save_every

nsave      = M0/save_every; % we store nsave elements
ind_save   = 0;             % counter for save


%==========================================================================
%% Setup info/settings and IDENTIFYING RESTRICTIONS
%==========================================================================
info=SetupInfo(nvar,m, nlag,horizons,@(x)chol(x));


%==========================================================================
%% Mappings
%==========================================================================
fo                 = @(x)f_h(x,info);
fo_inv             = @(x)f_h_inv(x,info);
fo_str2irfs        = @(x)StructuralToIRF(x,info);
fo_str2irfs_inv    = @(x)IRFToStructural(x,info);
gs_qr    = @(x)qr_unique(x);
%==========================================================================
%% write data in Rubio, Waggoner, and Zha (RES 2010)'s notation
%==========================================================================
% yt(t) A0 = xt(t) Aplus + constant + et(t) for t=1...,T;
% yt(t)    = xt(t) B     + ut(t)            for t=1...,T;
% x(t)     = [constant,yt(t-1), ... ,yt(t-nlag)];
% matrix notation yt = xt*B + ut;
% xt=[yt_{-1} ones(T,1)];

nobs0=8;
y0bar = mean(num(1:nobs0,:),1);% useful for GLP prior
info.y0bar = y0bar;

Y0 = data(1:8,var_id);  % save the first 8 obs as the initial conditions
Y = data(9:end,var_id);
[T,nvar] = size(Y);
tmpY = [Y0(end-nlag+1:end,:); Y];
X = zeros(T,nvar*nlag);
for ii=1:nlag
    X(:,(ii-1)*nvar+1:ii*nvar) = tmpY(nlag-ii+1:end-ii,:);
end
X = [ones(T,1) X];
%Z = X;




% find the optimal kappa values
[ml_opt,kappa] = get_OptKappa(Y0,Y,X,nlag,[.04,.0016],'redu',idx_ns);

sig2 = get_resid_var(Y0,Y);
prior_stru = prior_ACP_stru(nvar,nlag,kappa,sig2,idx_ns);
prior_redu = prior_ACP_redu(nvar,nlag,kappa,sig2,idx_ns);



%% posterior for reduced-form parameters
%[T,n] = size(Y);
% tmpY = [Y0(end-nlag+1:end,:); Y];
% Z = zeros(T,nvar*nlag);
% for ii=1:nlag
%     Z(:,(ii-1)*nvar+1:ii*nvar) = tmpY(nlag-ii+1:end-ii,:);
% end
% Z = [ones(T,1) Z];
k_beta = nvar^2*nlag+nvar;
k_alp = nvar*(nvar-1)/2;
%Beta = zeros(nsim,k_beta);
%Alp = zeros(nsim,k_alp);
%Sig = zeros(nsim,nvar);
count_alp = 0;
nsim=1000;

Thetai = cell([nvar,nsim]);
Sigi=nan(nvar,nsim);
% Create the cell array
Thetais = cell([nvar,1]);


for ii = 1:nvar
    yi = Y(:,ii);
    ki = nvar*nlag+ii;
    mi = [prior_redu.beta0(:,ii);prior_redu.alp0(count_alp+1:count_alp+ii-1)];
    Vi = sparse(1:ki,1:ki,[prior_redu.Vbeta(:,ii);prior_redu.Valp(count_alp+1:count_alp+ii-1)]);
    nui = prior_redu.nu(ii);
    Si = prior_redu.S(ii);
    Xi = [X -Y(:,1:ii-1)];
    % compute the parameters of the posterior distribution
    iVi = Vi\speye(ki);
    Kthetai = iVi + Xi'*Xi;
    CKthetai = chol(Kthetai,'lower');
    thetai_tilde = (CKthetai')\(CKthetai\(iVi*mi + Xi'*yi));
    Si_tilde = Si + (yi'*yi + mi'*iVi*mi - thetai_tilde'*Kthetai*thetai_tilde)/2;
    
    posterior_redu.nnu{ii}=prior_redu.nu(ii) + (T+ ki)/2;
    posterior_redu.Si_tilde{ii}=Si_tilde;
    posterior_redu.thetai_tilde{ii}=thetai_tilde;
    posterior_redu.Kthetai{ii} = Kthetai;
    posterior_redu.ki{ii} = size(Kthetai,1);
    
    posterior_redu.cholinvKthetai{ii} = chol(posterior_redu.Kthetai{ii}\eye(posterior_redu.ki{ii}),'lower');
    
    
    count_alp = count_alp + ii -1;
    
    % sample sig and theta
    
    
    
    for s=1:nsim
        Sigi_tmp=1./gamrnd(nui+T/2,1./Si_tilde,1,1);
        Sigi(ii,s) = Sigi_tmp;
        U = randn(1,ki).*repmat(sqrt(Sigi_tmp),1,ki);
        
        Thetai{ii,s} = repmat(thetai_tilde',1,1) + U/CKthetai;
    end
    
end

Sigmadraw_all=zeros(nvar,nvar,nsim);
Bdraw_all=zeros(m,nvar,nsim);
for s=1:nsim
    % Extract Thetai{ii,s} for each ii
    for ii = 1:nvar
        
        Thetais{ii} =  Thetai{ii,s};
    end
    
    [Bdraw_all_tmp,Sigmadraw_all_tmp] =   tthetassigma2TOBSigma(Thetais,Sigi(:,s),info);
    
    Bdraw_all(:,:,s) =Bdraw_all_tmp;
    Sigmadraw_all(:,:,s)=Sigmadraw_all_tmp;
    
end
mean_Sigmadraw_all = mean(Sigmadraw_all,3);
mean_Bdraw_all     = mean(Bdraw_all,3);



posterior_redu.Rmean    =  cell([nvar,1]);
posterior_redu.Rvariance = cell([nvar,1]);
for ii=1:nvar
    ki = nvar*nlag+ii;
    posterior_redu.Rmean{ii} = zeros(1,2*posterior_redu.nnu{ii});
    
end



%% useful definitions
% definitios used to store orthogonal-reduced-form draws, volume elements, and unnormalized weights
Bdraws         = cell([nsave,1]); % reduced-form lag parameters
Sigmadraws     = cell([nsave,1]); % reduced-form covariance matrices
Qdraws         = cell([nsave,1]); % orthogonal matrices


% definitions related to IRFs and stability of the coefficients; based on page 12 of Rubio, Waggoner, and Zha (RES 2010)
J      = [e;repmat(zeros(nvar),nlag-1,1)];
A      = cell(nlag,1);
extraF = repmat(zeros(nvar),1,nlag-1);
bigF      = zeros(nlag*nvar,nlag*nvar);
for l=1:nlag-1
    bigF((l-1)*nvar+1:l*nvar,nvar+1:nlag*nvar)=[repmat(zeros(nvar),1,l-1) e repmat(zeros(nvar),1,nlag-(l+1))];
end
info.extraF = extraF;
info.J = J;
info.A = A;
info.bigF = bigF;


% definition to facilitate the draws from B|Sigma
hh              = info.h;



BHat            = mean_Bdraw_all;
SigmaHat           =mean_Sigmadraw_all;
SigmaHat           = (SigmaHat'+SigmaHat)*0.5;




[ttheta_old,ssigma_old] = BSigmaTOtthetassigma2(BHat,SigmaHat);



% check inverse
[BHat2,SigmaHat2] = tthetassigma2TOBSigma(ttheta_old,ssigma_old,info);

% check inverse on the opposite direction
[ttheta_old2,ssigma_old2] = BSigmaTOtthetassigma2(BHat2,SigmaHat2);

% Usage:
identical = compare_cell_arrays(ttheta_old,ttheta_old2,1e-8);
if identical
    disp('Cell arrays are identical');
else
    disp('Cell arrays are NOT identical');
end


Sigmadraw=SigmaHat;
SigmaHatInv = inv(SigmaHat);
% % Check if SigmaHatInv is positive definite
% % Compute the Cholesky decomposition of SigmaHatInv
L = chol(SigmaHatInv, 'lower');
% % L*L' = SigmaHatInv, so R = L
R_base = chol( diag(ssigma_old)\eye(nvar), 'lower');
R_old = [R_base, zeros(nvar, T - nvar)];


init = 0;

[Q0,~]=DrawQ(info.nvar);
L0init = L*Q0;




info.nshocks=8;



% shock 1= demand
nr1=14;
info.Ss{1,1}= zeros(nr1,nvar);
info.Ss{1,1}(1,1)=1; % GDP +
info.Ss{1,1}(2,10)=1; % GDP DEF +
info.Ss{1,1}(3,11)=1; % PCE INDEX +
info.Ss{1,1}(4,12)=1; % PCEXFE INDEX +
info.Ss{1,1}(5,13)=1; % CPI +
info.Ss{1,1}(6,14)=1; % CPIXFE +
info.Ss{1,1}(7,19)=-1; % UNEMPLOYMENT -
info.Ss{1,1}(8,20)=1; % INDPRO +
info.Ss{1,1}(9,21)=1;  % CAP UTIL +
info.Ss{1,1}(10,25)=1; % FED FUNDS +
info.Ss{1,1}(11,26)=1; % 3-month t-bill +
info.Ss{1,1}(12,30)=1; % PRIME RATE +
info.Ss{1,1}(13,4)=-1;info.Ss{1,1}(13,1)=+1; % NON RESIDENTIAL minus GDP
info.Ss{1,1}(14,7)=-1;info.Ss{1,1}(14,1)=+1; % GOV SPENDING minus GDP



% shock 2=investment
nr2=15;
info.Ss{2,1}= zeros(nr2,nvar);
info.Ss{2,1}(1,1)=1;    % GDP +
info.Ss{2,1}(2,10)=1;   % GDP DEF +
info.Ss{2,1}(3,11)=1;   % PCE index +
info.Ss{2,1}(4,12)=1;   % PCEXFE index +
info.Ss{2,1}(5,13)=1;   % CPI +
info.Ss{2,1}(6,14)=1;   % CPIXFE +
info.Ss{2,1}(7,19)=-1;  % U -
info.Ss{2,1}(8,20)=1;   % INDPRO +
info.Ss{2,1}(9,21)=1;   % CAP UTIL +
info.Ss{2,1}(10,25)=1;  % FED FUNDS +
info.Ss{2,1}(11,26)=1;  % 3-MONTH T-BILL +
info.Ss{2,1}(12,30)=1;  % PRIME RATE +
info.Ss{2,1}(13,34)=-1;  % S&P 500 -
info.Ss{2,1}(14,4)=1;info.Ss{2,1}(14,1)=-1;
info.Ss{2,1}(15,7)=-1;info.Ss{2,1}(15,1)=+1;

% % shock 3=financial
nr3=15;
info.Ss{3,1}= zeros(nr3,nvar);
info.Ss{3,1}(1,1)=1;    % GDP +
info.Ss{3,1}(2,10)=1;   % GDP DEF +
info.Ss{3,1}(3,11)=1;   % PCE index +
info.Ss{3,1}(4,12)=1;   % PCEXFE index +
info.Ss{3,1}(5,13)=1;   % CPI +
info.Ss{3,1}(6,14)=1;   % CPIXFE +
info.Ss{3,1}(7,19)=-1;  % U -
info.Ss{3,1}(8,20)=1;   % INDPRO +
info.Ss{3,1}(9,21)=1;   % CAP UTIL +
info.Ss{3,1}(10,25)=1;  % FED FUNDS +
info.Ss{3,1}(11,26)=1;  % 3-MONTH T-BILL +
info.Ss{3,1}(12,30)=1;  % PRIME RATE +
info.Ss{3,1}(13,34)=1;  % S&P 500 +
info.Ss{3,1}(14,4)=1; info.Ss{3,1}(14,1)=-1;
info.Ss{3,1}(15,7)=-1;info.Ss{3,1}(15,1)=+1;


% shock 4=monetary policy
nr4=19;
info.Ss{4,1}= zeros(nr4,nvar);
info.Ss{4,1}(1,1)=-1;    % GDP -
info.Ss{4,1}(2,10)=-1;   % GDP DEF -
info.Ss{4,1}(3,11)=-1;   % PCE index -
info.Ss{4,1}(4,12)=-1;   % PCEXFE index -
info.Ss{4,1}(5,13)=-1;   % CPI -
info.Ss{4,1}(6,14)=-1;   % CPIXFE -
info.Ss{4,1}(7,18)=-1;  % EMPL -
info.Ss{4,1}(8,19)= 1;  % U +
info.Ss{4,1}(9,20)=-1;   % INDPRO -
info.Ss{4,1}(10,21)=-1;   % CAP UTIL -
info.Ss{4,1}(11,25)=1;  % FED FUNDS +
info.Ss{4,1}(12,26)=1;  % 3-MONTH T-BILL +
info.Ss{4,1}(13,27)=1;  % 2-YEAR T-NOTE +
info.Ss{4,1}(14,28)=1;  % 5-YEAR T-NOTE +
info.Ss{4,1}(15,29)=1;  % 10-YEAR T-NOTE +
info.Ss{4,1}(16,30)=1;  % PRIME RATE +
info.Ss{4,1}(17,31)=1;  % AAA CORPORTATe +
info.Ss{4,1}(18,32)=1;  % BAA CORPORTATe +
info.Ss{4,1}(19,34)=-1;  % S&P 500 +

%
%
% % shock 5= government
nr5=14;
info.Ss{5,1}= zeros(nr5,nvar);
info.Ss{5,1}(1,1)= 1; % GDP +
info.Ss{5,1}(2,7)= 1; % GOV SPENDING +
info.Ss{5,1}(3,8)= -1; % FED BUDGET +
info.Ss{5,1}(4,9)= 1; % FED TAX RECEIPTS +
info.Ss{5,1}(5,10)=1;   % GDP DEF +
info.Ss{5,1}(6,11)=1;   % PCE index +
info.Ss{5,1}(7,12)=1;   % PCEXFE index +
info.Ss{5,1}(8,13)=1;   % CPI +
info.Ss{5,1}(9,14)=1;   % CPIXFE +
info.Ss{5,1}(10,19)=-1; % UNEMPLOYMENT -
info.Ss{5,1}(11,25)=1; % FFR +
info.Ss{5,1}(12,26)=1; % 3-MONTH T-BILL +
info.Ss{5,1}(13,30)=1; % PRIME +
info.Ss{5,1}(14,7)=1;info.Ss{5,1}(14,1)=-1;
government = [1,NaN,NaN,NaN,NaN,NaN,1,-1,1,1,1,1,1,1,NaN,NaN,NaN,NaN,NaN,...
    NaN,NaN,NaN,NaN,NaN,1,1,NaN,NaN,NaN,1,NaN,NaN,NaN,NaN,NaN]';
%
%
% shock 6= technology
nr6=12;
info.Ss{6,1}= zeros(nr6,nvar);
info.Ss{6,1}(1,1)= 1; % GDP +
info.Ss{6,1}(2,2)= 1; % CONS EXPENDITURE +
info.Ss{6,1}(3,4)= 1; % NON RESIDENTIAL +
info.Ss{6,1}(4,10)=-1;   % GDP DEF -
info.Ss{6,1}(5,11)=-1;   % PCE index -
info.Ss{6,1}(6,12)=-1;   % PCEXFE index -
info.Ss{6,1}(7,13)=-1;   % CPI -
info.Ss{6,1}(8,14)=-1;   % CPIXFE -
info.Ss{6,1}(9,15)= 1;   % HOURLY WAGE +
info.Ss{6,1}(10,16)= 1;   % LABOR PROD +
info.Ss{6,1}(11,17)= 1;   % UTILIZATION ADJUSTED +
info.Ss{6,1}(12,19)=-1; % UNEMPLOYMENT -
technology = [1,1,NaN,1,NaN,NaN,NaN,NaN,NaN,-1,-1,-1,-1,-1,1,1,1,NaN,NaN,...
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]';
%
% % shock 7= labor supply
nr7=8;
info.Ss{7,1}= zeros(nr7,nvar);
info.Ss{7,1}(1,1)= 1; % GDP +
info.Ss{7,1}(2,10)=-1;   % GDP DEF -
info.Ss{7,1}(3,11)=-1;   % PCE index -
info.Ss{7,1}(4,12)=-1;   % PCEXFE index -
info.Ss{7,1}(5,13)=-1;   % CPI -
info.Ss{7,1}(6,14)=-1;   % CPIXFE -
info.Ss{7,1}(7,15)= -1;   % HOURLY WAGE -
info.Ss{7,1}(8,19)= 1; % UNEMPLOYMENT +
labor = [1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,-1,-1,-1,-1,-1,-1,NaN,NaN,NaN,1,...
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]';
sum(isnan(labor)~=1)
% shock 8= wage bargaining
nr8=8;
info.Ss{8,1}= zeros(nr8,nvar);
info.Ss{8,1}(1,1)= 1; % GDP +
info.Ss{8,1}(2,10)=-1;   % GDP DEF -
info.Ss{8,1}(3,11)=-1;   % PCE index -
info.Ss{8,1}(4,12)=-1;   % PCEXFE index -
info.Ss{8,1}(5,13)=-1;   % CPI -
info.Ss{8,1}(6,14)=-1;   % CPIXFE -
info.Ss{8,1}(7,15)= -1;   % HOURLY WAGE -
info.Ss{8,1}(8,19)= -1; % UNEMPLOYMENT -
wage = [1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,-1,-1,-1,-1,-1,-1,NaN,NaN,NaN,-1,...
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]';






% shock 9= oil
nr9=13;
info.Ss{9,1}= zeros(nr9,nvar);
info.Ss{9,1}(1,1)= 1; % GDP +
info.Ss{9,1}(2,2)= 1; % CONS EXPENDITURE +
info.Ss{9,1}(3,4)= 1; % NON RESIDENTIAL +
info.Ss{9,1}(4,10)=-1;   % GDP DEF -
info.Ss{9,1}(5,11)=-1;   % PCE index -
info.Ss{9,1}(6,12)=-1;   % PCEXFE index -
info.Ss{9,1}(7,13)=-1;   % CPI -
info.Ss{9,1}(8,14)=-1;   % CPIXFE -
info.Ss{9,1}(9,15)= 1;   % HOURLY WAGE +
info.Ss{9,1}(10,16)= 1;   % LABOR PROD +
info.Ss{9,1}(11,17)= 1;   % UTILIZATION ADJUSTED +
info.Ss{9,1}(12,19)= 1; % UNEMPLOYMENT
info.Ss{9,1}(13,35)= -1; % oil price


% shock 10= sentiment

nr1=11;
info.Ss{10,1}= zeros(nr1,nvar);
info.Ss{10,1}(1,1)=1; % GDP +
info.Ss{10,1}(2,2)=1; % PCE +
info.Ss{10,1}(3,3)=1; % residential +
info.Ss{10,1}(4,4)=1; % non-residential +
info.Ss{10,1}(5,10)=1; % GDP DEF +
info.Ss{10,1}(6,11)=1; % PCE INDEX +
info.Ss{10,1}(7,12)=1; % PCEXFE INDEX +
info.Ss{10,1}(8,13)=1; % CPI +
info.Ss{10,1}(9,14)=1; % CPIXFE +
info.Ss{10,1}(10,19)=1; % UNEMPLOYMENT -
info.Ss{10,1}(11,34)=1; % s&p +



Bdraw=BHat;
z_old_R=cell([nvar,1]);
for i=1:nvar
    z_old_R{i} = zeros(1,2*posterior_redu.nnu{i})';
    z_old_R{i}(1,1) = R_old(i,i);
end




% =====================================
% New initialization
% =====================================
disp('initilization of Q ...');
Qdraw = draw_Q_conditional_on_sign(Sigmadraw, info);
disp('initilization of Q ... done ...');


% Create a sample upper triangular matrix R
n = size(Qdraw, 1);
Rdraw_check = triu(rand(n, n));  % Random upper triangular matrix

% Compute X
XX_check = Qdraw * Rdraw_check;


%% initialize counters to track the state of the computations
counter = 1;
count   = 0;

record = 1;






switch label_R
    
    case 'nolabel'
        
        function_restrictions = @restriction_none;
        
    case {'cmy'}
        
        function_restrictions = @restriction_cmy;
        
    otherwise
        disp('choose a valid specification')
        return
end


z_old_X =vec(XX_check);
z_old_B =vec(Bdraw);


% check restrictions

[Qdraw,~]=gs_qr(XX_check);

function_restrictions_i = @(ag0, ag1, ag2, ag3, ag4, ag5) function_restrictions(ag0, ag1, ag2, ag3, ag4, ag5);

S_prop = function_restrictions_i(Bdraw, Sigmadraw, Qdraw,fo_inv,fo_str2irfs,info);

Bold = Bdraw; %reduced form
Sigold = Sigmadraw;
Ainvold = chol(Sigold,'lower') /diag(sqrt(ssigma_old));
Aold = eye(size(Sigold))/Ainvold;
Bstrold = Bold*Aold';

S_prop = function_restrictions(Bold, Sigold, Qdraw,fo_inv,fo_str2irfs,info);

clear Bdraw Sigmadraw




%posterior_Sigmadraws_accepted = nan(M0,1);
tstart = tic;
while record<=M0
    
    %record
    
    %% setting function restrictions
    function_restrictions_i = @(ag0, ag1, ag2, ag3, ag4, ag5) function_restrictions(ag0, ag1, ag2, ag3, ag4, ag5);
    
    
    %% Draw Q
    
    slice.scale_z    = 1;
    slice.mean       = zeros(info.nvar*info.nvar,1);
    slice.chol_cov_z = chol(eye(info.nvar*info.nvar),'lower');
    slice.chol_cov_z = 1;
    slice.fcn_lik    = @(z_prop) loglike_Q(z_prop,function_restrictions_i,gs_qr,fo_inv,fo_str2irfs,Bold,Sigold,info);
    slice.nobs       = info.nvar*info.nvar;
    
    % --- actual slice sampling
    % z_old = randn(info.nvar*info.nvar,1);
    lik_old = slice.fcn_lik(z_old_X);
    [z_old_X, ~, n_try] = slice_sampling_v02(slice, z_old_X, lik_old);
    [Qdraw,Rdraw] =  gs_qr(reshape(z_old_X,info.nvar,info.nvar));
    
    %% Draw sigma^2=(sigma_1^2,\dots,\sigma_n^2)
    % draw equation by equation
    %
    for ii=1:nvar
        slice.scale_z    = 1;
        slice.mean       = vec(posterior_redu.Rmean{ii});
        
        nnuTilde_i = 2*posterior_redu.nnu{ii};
        Stilde_ii =posterior_redu.Si_tilde{ii} + 0.5*(ttheta_old{ii}-posterior_redu.thetai_tilde{ii})'*posterior_redu.Kthetai{ii}*(ttheta_old{ii}-posterior_redu.thetai_tilde{ii});
        posterior_redu.Rvariance{ii} = 1/(2*Stilde_ii);% depends on theta
        chol_lower_Rvariance_common = chol(posterior_redu.Rvariance{ii}, 'lower');
        slice.chol_cov_z_common = chol_lower_Rvariance_common;slice.nobs1 = 1;
        slice.nobs2 = nnuTilde_i;
        slice.fcn_lik    = @(z_prop) loglike_sigma2i_mc_BSig_v02(z_prop,function_restrictions_i,fo_inv,fo_str2irfs,ttheta_old,ssigma_old,Qdraw,info,nnuTilde_i,ii,Bold,Sigold,Aold,Bstrold);
        slice.nobs       = slice.nobs1*slice.nobs2;
        
        % --- actual slice sampling        
        lik_old = slice.fcn_lik(z_old_R{ii});
        
        [z_old_R_i, ~, n_try, Bold, Sigold, Aold,Bstrold] = slice_sampling_v02_Sigma_BSig(slice, z_old_R{ii}, lik_old);
        z_old_R{ii} = z_old_R_i;
        
        ssigma_old(ii,1)=1/(z_old_R_i'*z_old_R_i);

    end

    % %% Draw B
    % % Draw \theta
    % % draw equation by equation
    for ii=1:nvar
        cholLSigmadraw =  chol(ssigma_old(ii,1));

        slice.scale_z    = 1;
        slice.mean       = vec(posterior_redu.thetai_tilde{ii});
        
        slice.chol_cov_z1 = cholLSigmadraw;
        slice.chol_cov_z2 = posterior_redu.cholinvKthetai{ii};
        slice.nobs1 = 1;
        slice.nobs2 = posterior_redu.ki{ii};
        
        slice.fcn_lik    = @(z_prop) loglike_tthetai_BSig(z_prop,function_restrictions_i,fo_inv,fo_str2irfs,ttheta_old,ssigma_old,Qdraw,info,ki,ii,Bold,Sigold,Aold,Bstrold);
        
        slice.nobs       = slice.nobs1*slice.nobs2;
        
        % --- actual slice sampling
        
        lik_old = slice.fcn_lik(ttheta_old{ii});
        
        %[z_old_B, ~, n_try] = slice_sampling_v02(slice, z_old_B, lik_old);
        [ttheta_old_i, ~, n_try,Bold,Sigold,Aold,Bstrold] = slice_sampling_v02_kron_BSig(slice, ttheta_old{ii}, lik_old);
        
        ttheta_old{ii} = ttheta_old_i;
        
    end
    
    % [Bdraw,Sigmadraw] = tthetassigma2TOBSigma(ttheta_old,ssigma_old,info);

    %% check if sign restrictions hold
    
    %   L0 =reshape(irfpara(1:nvar*nvar,:),nvar,nvar);
    
    %    BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];

    count=count+1;    
    if rem(record, save_every) == 0
        % --- store
        % store orthogonal reduced-form draws
        ind_save = ind_save + 1;
        Bdraws{ind_save,1}     = Bold;
        Sigmadraws{ind_save,1} = Sigold;
        Qdraws{ind_save,1} = Qdraw;
    end

    record=record+1;

    if counter==iter_show
        
        display(['Number of draws = ',num2str(record)])
        display(['Remaining draws = ',num2str(M0-(record))])
        counter =0;
        
    end
    counter = counter + 1;
    
end

telapsed = toc(tstart);

disp(telapsed)

%% store draws
L    = zeros(horizon+1,nvar,info.nshocks,ind_save);
%cumL    = zeros(horizon+1,nvar,1:,ind_save);


for s=1:ind_save
    
    
    
    Bdraw =     Bdraws{s,1} ;
    Sigmadraw = Sigmadraws{s,1} ;
    Qdraw=Qdraws{s,1};
    
    BSigmaQ = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];
    
    
    
    structpara = f_h_inv(BSigmaQ,info);
    
    
    LIRF = IRF_horizons(structpara, nvar, nlag, m, nex, 0:horizon);
    
    
    
    
    
    for h=0:horizon
        L(h+1,:,1:info.nshocks,s) =  LIRF(1+h*nvar:(h+1)*nvar,1:info.nshocks);
        
        % for i=1:nvar
        %     cumL(h+1,1:4,i,s)   = sum(L(1:h+1,1:4,i,s),1);
        % end
        
    end
    
    
    
end

% stack = cat(3, Sigmadraws{:});        % 10×3×3 array
% SigmaMean = mean(stack, 3);
% 
% stack = cat(3, Bdraws{:});
% BMean = mean(stack, 3);


cd results

if ~exist('matfiles', 'dir')
    mkdir('matfiles')
end


savefile= ['matfiles/rgibbs_',label_R,'asymmetric_prior__ndraws',num2str(M0),'n_shocks',num2str(info.nshocks),'.mat'];

save(savefile,'L','horizon','telapsed','-v7.3');
cd ..
