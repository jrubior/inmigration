%==========================================================================
%% housekeeping
%==========================================================================
clear variables;close all;userpath('clear');restoredefaultpath;clc;
tic;
rng('default'); % reinitialize the random number generator to its startup configuration
seed =0;
rng(seed,'twister');         % set seed
currdir=pwd;
addpath([currdir,'/data']); % set path to data functions
addpath([currdir,'/helpfunctions']); % set path to helper functions
addpath([currdir,'/helpfunctions/ChrisSimsOptimize']);
addpath([currdir,'/helpfunctions/subroutines']); % set path to helper functions
cd(currdir)


%==========================================================================
%% load the data
%==========================================================================
data = readmatrix('data35_1975.xlsx','Range','B9:AJ185');
        % 35 variables + 8 shocks
        % shocks: Demand, Investment, Financial, Monetary, Govt spending
        %   Technology, Labor supply, Wage bargaining
    var_id = 1:35;
    idx_ns = 1:35; % index for variables in levels   

num = data(:,var_id);

%=========================================================================
%% model setup
%==========================================================================
nlag      = 5;                     % number of lags
nvar      = size(num,2);           % number of endogenous variables
nex       = 1;                     % 12 because of seasonal dummies
m         = nvar*nlag + nex;       % number of exogenous variables
horizon   = 35;                    % maximum horizon for IRFs
horizons  = [0,1,2,3,4,inf];       % horizons upon which sign and zero restrictions can be imposed
NS        = 1 + numel(horizons);   % number of objects in F(THETA) to which we impose sign and zero restrictions: F(THETA)=[A_0;L_{0};L_{1};L_{2};L_{3};L_{4};L_{inf}]
e         = eye(nvar);             % create identity matrix
M0  =5e2;
% M0 = 1e4;
conjugate  = 'none';
iter_show  = 1e2;
fixed_rf   = 0;
prior_only = 0;
label_R = 'cmy';
prior_type = 'minnesota';
if fixed_rf==1
    draw_Sigma=0;
    draw_B=0;
else
    draw_Sigma=1;
    draw_B=1;
end
reoptimizeMH=0;

save_every   = 1000;%50;%50;%50;             % save every save_every

nsave      = M0/save_every; % we store nsave elements
ind_save   = 0;             % counter for save

elliptical_Q = 1; %0 when using Chan, Matthes, Yu, 1 when using elliptical
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


Ynd = Y;
Xnd = X;


%% prior for reduced-form parameters



switch prior_type
    case 'flat'


    %% prior
    nnuBar             = 0;
    OomegaBarInverse   = zeros(m);
    mmuBar             = zeros(m,nvar);  % Psi in the paper
    PpsiBar            = zeros(nvar);    % Phi in the paper

    case 'minnesota'

        % Giannone, Lenza, Primiceri (2015) hyperpriors modes
        mode.llambda = .2;

        % Giannone, Lenza, Primiceri (2015) hyperpriors std deviations
        sd.llambda = .4;

        % Giannone, Lenza, Primiceri (2015) scale and shape of the IG on psi/(d-n-1)
        scalePSI   = 0.02^2;    


         % transform mode and std into Gamma(k,theta) coefficients
        info.priorcoef.llambda= GammaCoef(mode.llambda,sd.llambda,0);
        info.priorcoef.alpha.PSI=scalePSI;
        info.priorcoef.beta.PSI=scalePSI;
        info.pos = []; % stationary variables


        % bounds for maximization
        MIN.llambda = 0.0001;
        MAX.llambda = 5;



        % initial guess for hyper-parameters
        llambda  = 0.5;
        info.Vc  = 10^7;

        % initial guess for psi (phi in the paper)
        % residual variance of AR(1) for each variable
        SS=zeros(nvar,1);
        for i=1:nvar
            yols=Y(2:end,i);
            xols=[ones(T-1,1),Y(1:end-1,i)];
            bbetaols = (xols'*xols)\(xols'*yols);
            eols = yols-xols*bbetaols;
            SS(i)=eols'*eols/(size(yols,1)-size(xols,2));
        end
        MIN.ppsi  = SS./100;
        MAX.ppsi  = SS.*100;
        ppsi      = SS;
        inllambda = -log((MAX.llambda-llambda)./(llambda-MIN.llambda));
        inppsi    = -log((MAX.ppsi-ppsi)./(ppsi-MIN.ppsi));

        x0 = [inllambda;inppsi];


        f0 = negative_loglike_obj_mp(x0,info,Xnd,Ynd,MIN,MAX);

        crit                              = 1e-16;
        nit                               = 1000;
        H0                                = eye(size(x0,1))*10;


        if reoptimizeMH==1

            [fh,xh,~,~,itct,~,~]              = csminwel('negative_loglike_obj_mp',x0,H0,[],crit,nit,info,Xnd,Ynd,MIN,MAX);

            save('xh.mat','xh','fh');
        else
            load('xh.mat','xh','fh');
        end

        llambda0_xh = MIN.llambda+(MAX.llambda-MIN.llambda)./(1+exp(-xh(1,1)));
        ppsi0_xh    = MIN.ppsi+(MAX.ppsi-MIN.ppsi)./(1+exp(-xh(2:info.nvar+1,1)));%
  
        llambda0 =  llambda0_xh;
        ppsi0    =  ppsi0_xh;


       % we re-set the random seed after using csminwel to guarantee the
            % posterior draws are identical across operating systems.
        rng(seed,'twister');
   


        % set position where variables are in first differences

        % inverse Wishart prior parameters
        nnuBar              = info.nvar + 2;
        PpsiBar             = diag(ppsi0);
        oomegaBar=zeros(info.m,1);
        oomegaBar(1:info.nex,1)=info.Vc;
        for ell=1:info.nlag
            oomegaBar(info.nex+1+(ell-1)*info.nvar:info.nex+ell*info.nvar,1) =  (llambda0^2)*(nnuBar-info.nvar-1)./((ell^2)*ppsi0);
        end
        OomegaBar=diag(oomegaBar);
        OomegaBarInverse=diag(1./oomegaBar);
        mmuBar = zeros(info.m,info.nvar);
        diagmmuBar=ones(info.nvar,1);
        diagmmuBar(info.pos)=0;   % Set to zero the prior mean on the first own lag for variables selected in the vector pos
    
        mmuBar(info.nex+1:info.nvar+info.nex,:)=diag(diagmmuBar);


        Y = Ynd;
        X = Xnd;
        T = size(Y,1);


end




%% posterior for reduced-form parameters
nnuTilde            = T +nnuBar;
OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(m);
OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
mmuTilde            = (X'*X  + OomegaBarInverse)\(X'*Y + OomegaBarInverse*mmuBar);
PpsiTilde           = Y'*Y + PpsiBar + mmuBar'*OomegaBarInverse*mmuBar - mmuTilde'*OomegaTildeInverse*mmuTilde;
PpsiTilde           = (PpsiTilde'+PpsiTilde)*0.5;

%% Parameters to fix reduced-form
BHat            = mmuTilde;%(X'*X)\(X'*Y);
SigmaHat           =PpsiTilde/(nnuTilde+size(PpsiTilde,1)+1);% (Y-X*BHat)'*(Y-X*BHat)/size(Y,1);
SigmaHat           = (SigmaHat'+SigmaHat)*0.5;


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
cholOomegaTilde = hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below

chol_lower_OomegaTilde = chol(OomegaTilde, 'lower');
chol_lower_OomegaTilde_inv = eye(size(OomegaTilde))/chol_lower_OomegaTilde;

if prior_only==1
    cholOomegaBar   = hh(OomegaBar)'; % this matrix is used to draw B|Sigma below
end


Rmean = nan(nvar,nnuTilde);
for i=1:nnuTilde
    Rmean(:,i) = zeros(nvar,1);
end
Rvariance = zeros(nvar*nnuTilde,nvar*nnuTilde);
for i=1:nnuTilde
    Rvariance((i-1)*nvar+1:i*nvar,(i-1)*nvar+1:i*nvar) = inv(PpsiTilde);
end
chol_lower_Rvariance = chol(Rvariance,'lower');
chol_lower_Rvariance_common = chol(inv(PpsiTilde), 'lower');									

% initialize Gibbs
Sigmadraw=SigmaHat;
SigmaHatInv = inv(SigmaHat);
% Check if SigmaHatInv is positive definite
% Compute the Cholesky decomposition of SigmaHatInv
L = chol(SigmaHatInv, 'lower');
% L*L' = SigmaHatInv, so R = L
R_base = L;
R = [R_base, zeros(nvar, nnuTilde - nvar)];

init = 0;

[Q0,~]=DrawQ(info.nvar);
L0init = L*Q0;




info.nshocks=10;



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

% index=find(info.Ss{1,1}*L0init(:,1)<0)
% varindex=find((sum(info.Ss{1,1}(index,:),1)~=0)==1)
% L0init(varindex,1)=-L0init(varindex,1);
% indexnew=find(info.Rs{1,1}*L0init(:,1)<0)


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



% L0initNEW=L0init;
% for j=1:info.nshocks
% index=find(info.Ss{j,1}*L0init(:,j)<0)
% varindex=find((sum(info.Ss{j,1}(index,:),1)~=0)==1)
% L0initNEW(varindex,j)=-L0init(varindex,j);
% ndexnew=find(info.Ss{j,1}*L0initNEW(:,j)<0)
% end
% 
% Sigmadraw= L0initNEW*L0initNEW';
% Sigmadraw = (Sigmadraw+Sigmadraw')/2;
% Qdraw = chol(Sigmadraw,'lower')\L0initNEW;
% 
% 
% SigmadrawInv = inv(Sigmadraw);
% % Check if SigmaHatInv is positive definite
% % Compute the Cholesky decomposition of SigmaHatInv
% L = chol(SigmadrawInv, 'lower');
% % L*L' = SigmaHatInv, so R = L
% R_base = L;
% R = [R_base, zeros(nvar, nnuTilde - nvar)];




Bdraw=BHat;
z_old_R = vec(R);
  

%SigmaCMY     = load('SigmaCMY.mat','Sigtilde','Qdraw');
%Sigmadraw = SigmaCMY.Sigtilde;
%QdrawCMY = SigmaCMY.Qdraw;

% ========================================================================
% OLD initialization (gets slower when the number of shocks increases)
% ========================================================================

% tic
% 
% while init==0
% 
% 
% %     z_old_R = chol(Rvariance,'lower')*randn(info.nvar*nnuTilde,1);
% % 
% %     R_old = reshape(z_old_R,info.nvar,nnuTilde);
% %     Sigmadraw=inv(R_old*R_old');
%  % Sigmadraw     = iwishrnd(PpsiTilde,nnuTilde);
% 
% 
% 
%     Xold =randn(info.nvar);
%     [Qold,Rold]  = qr(Xold);
%     for ii=1:nvar
%         if Rold(ii,ii)<0
%             Qold(:,ii)=-Qold(:,ii);
%         end
%     end
% 
%    % Qold=QdrawCMY;
%     Rr=chol(Sigmadraw)'*Qold;
% 
% 
%     Tt=zeros(info.nshocks,nvar);
%     for j=1:info.nshocks
% 
%                 for i=1:nvar
%                     if all(info.Ss{j,1}*Rr(:,i)>0) 
%                         Tt(j,i)=1;
%                     elseif all(info.Ss{j,1}*Rr(:,i)<0) 
%                         Tt(j,i)=-1;
%                     end
%                 end
% 
% 
%         % switch j
%         %     case 1
%         %         for i=1:nvar
%         %             if all(info.Ss{j,1}*Rr(:,i)>0)
%         % 
%         %                 Tt(j,i)=1;
%         %             elseif all(info.Ss{j,1}*Rr(:,i)<0)
%         %                 Tt(j,i)=-1;
%         %             end
%         %         end
%         %     case {2,3,4,5,6,7,8}
%         % 
%         %         for i=1:nvar
%         %             if all(info.Ss{j,1}*Rr(:,i)>0) 
%         %                 Tt(j,i)=1;
%         %             elseif all(info.Ss{j,1}*Rr(:,i)<0) 
%         %                 Tt(j,i)=-1;
%         %             end
%         %         end
%         % end
% 
% 
%     end
% 
%     % checkA1 = any(sum(Tt ~= 0, 1) > 1);
%     % 
%     % if checkA1==1
%     %     disp('PROBLEM with Tt')
%     %     return;
%     % end
% 
% 
%     %if rank(Tt)<nshocks
% 
%     if any(sum(Tt == 0, 2)==size(Tt,2))
% 
%     else
% 
% 
%         Rstar = zeros(nvar,nvar);
%         S = cell(info.nshocks,1);
% 
%         isamples = zeros(info.nvar,1);
% 
%         % For each row j
%         for j = 1:info.nshocks
% 
%             S{j} = find(Tt(j,:) == 1 | Tt(j,:) == -1);
% 
% 
% 
%             % Get a random index from S_j
%             if numel(S{j})>1
%                 randomIndex = randsample(S{j},1);
%             else
%                 randomIndex = S{j};
%             end
% 
%             isamples(j,1) = randomIndex;
% 
% 
% 
%             if Tt(j,isamples(j))==1
%                 Rstar(:,j) = Rr(:,isamples(j,1));
%             elseif Tt(j,isamples(j))==-1
%                 Rstar(:,j) = -Rr(:,isamples(j,1));
%             end
% 
%         end
% 
%         % Loop through j = m+1, ..., n
%         %i_values=zeros(nvar-info.nshocks,1);
%         for j = info.nshocks+1:nvar
%             % Create S_j = {1,...,n} \ {i_1,...,i_j}
%             % We need to exclude i_1 through i_(j-1) from the set {1,...,n}
% 
% 
%             S_j = setdiff(1:nvar, isamples(1:j-1,1));
% 
% 
%             % Sample an element uniformly from S_j
%             if numel(S_j)>1
% 
% 
%                 random_index = randsample(S_j,1);
% 
%             else
%                 random_index = S_j;
%             end
%             sampled_element = random_index;
% 
% 
% 
%             % Store the sampled element as i_j
%             isamples(j,1) = sampled_element;
% 
%             if rand(1)>0.5
%                 Rstar(:,j) = Rr(:,isamples(j,1));
%             else
%                 Rstar(:,j) = -Rr(:,isamples(j,1));
%             end
% 
%         end
% 
%         Qdraw=(chol(Sigmadraw)')\Rstar;
% 
%         L0old=Rstar;
% 
% 
% 
% 
%         init = 1;
%     end
% end
% % Assuming Qdraw is already defined as an orthogonal matrix of size n×n
% init
% disp('time to find Q using CMY')
% toc

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

% % Verify using QR decomposition
% [Q_check, R_check] = qr(XX_check);
% keyboard
% % Adjust for potential sign differences
% for j = 1:n
%     if sign(Qdraw(1,j)) ~= sign(Q_check(1,j))  % Check first element of each column
%         Q_check(:,j) = -Q_check(:,j);      % Flip signs if different
%         R_check(j,:) = -R_check(j,:);
%     end
% end


%Bdraw=BHat;

%BSigmaQold = [vec(Bdraw);vec(SigmaHat);vec(Qold)];
%SigmaQold  = [vec(SigmaHat);vec(Qold)];



%% initialize counters to track the state of the computations
counter = 1;
count   = 0;

record = 1;






switch label_R

    case 'nolabel'

        function_restrictions = @restriction_none;

    case {'cmy','cmy10shocks'}

        %switch info.nshocks
            
         %   case 2
        %function_restrictions = @restriction_cmy_2shocks;
         %   case 3
         %      function_restrictions = @restriction_cmy_3shocks; 
        %            case 8
        function_restrictions = @restriction_cmy; 
       % end
        function_restrictions_Q = @restriction_cmy_withT;

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


%posterior_Sigmadraws_accepted = nan(M0,1);
tstart = tic;
while record<=M0
   % record
    if fixed_rf==1
        Sigmadraw     = SigmaHat;
        Bdraw         = BHat ;
    elseif prior_only==1
        %% Draw Sigma
        Sigmadraw     = iwishrnd(PpsiBar,nnuBar);
        cholLSigmadraw = hh(Sigmadraw)';
        Bdraw         = kron(cholLSigmadraw,cholOomegaBar)*randn(m*nvar,1) + reshape(mmuBar,nvar*m,1);
        Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);
    else


        % % if counter_Q==0
        %
        %      %counter_Q = counter_Q +1;
        %      %% Draw Sigma
        %      Sigmadraw     = iwishrnd(PpsiTilde,nnuTilde);
        %      cholSigmadraw = hh(Sigmadraw)';
        %      Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*nvar,1) + reshape(mmuTilde,nvar*m,1);
        %      Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);
        %  %end

    end


    %% setting function restrictions
    function_restrictions_i = @(ag0, ag1, ag2, ag3, ag4, ag5) function_restrictions(ag0, ag1, ag2, ag3, ag4, ag5);


    %% Draw Q

    if elliptical_Q==1
            
        slice.scale_z    = 1;
        slice.mean       = zeros(info.nvar*info.nvar,1);
        slice.chol_cov_z = chol(eye(info.nvar*info.nvar),'lower');
		slice.chol_cov_z = 1;

        slice.fcn_lik    = @(z_prop) loglike_Q(z_prop,function_restrictions_i,gs_qr,fo_inv,fo_str2irfs,Bdraw,Sigmadraw,info);
        slice.nobs       = info.nvar*info.nvar;

        % --- actual slice sampling
        % z_old = randn(info.nvar*info.nvar,1);
        lik_old = slice.fcn_lik(z_old_X);


        [z_old_X, ~, n_try] = slice_sampling_v02(slice, z_old_X, lik_old);



        [Qdraw,Rdraw] =  gs_qr(reshape(z_old_X,info.nvar,info.nvar));

    else

        S_prop=0;
        function_restrictions_Q_i = @(ag0, ag1, ag2, ag3, ag4, ag5) function_restrictions_Q(ag0, ag1, ag2, ag3, ag4, ag5);

        while S_prop==0


            [Q_prop,~]=DrawQ(info.nvar);

            [S_prop,Qdraw] = function_restrictions_Q_i(Bdraw, Sigmadraw, Q_prop,fo_inv,fo_str2irfs,info);

    
        end
      
          %S_prop = function_restrictions_i(Bdraw, Sigmadraw, Qdraw,fo_inv,fo_str2irfs,info);


    end

  

    %% Draw Sigma
    if draw_Sigma==1
        slice.scale_z    = 1;



        slice.mean       = vec(Rmean);
        %slice.chol_cov_z = chol(Rvariance,'lower');
        %slice.chol_cov_z = chol_lower_Rvariance; %
		slice.chol_cov_z_common = chol_lower_Rvariance_common;slice.nobs1 = info.nvar;
        slice.nobs2 = nnuTilde;						
        %slice.fcn_lik    = @(z_prop) loglike_Sigma(z_prop,function_restrictions_i,fo_inv,fo_str2irfs,Bdraw,Qdraw,info,nnuTilde,mmuTilde,OomegaTilde);
        slice.fcn_lik    = @(z_prop) loglike_Sigma_mc_v02(z_prop,function_restrictions_i,fo_inv,fo_str2irfs,Bdraw,Qdraw,info,nnuTilde,mmuTilde,chol_lower_OomegaTilde,chol_lower_OomegaTilde_inv);

        slice.nobs       = info.nvar*nnuTilde;

        % --- actual slice sampling

        lik_old = slice.fcn_lik(z_old_R);

        [z_old_R, ~, n_try] = slice_sampling_v02_Sigma(slice, z_old_R, lik_old);

        R_old = reshape(z_old_R,info.nvar,nnuTilde);
        Sigmadraw=inv(R_old*R_old');

       % posterior_Sigmadraws_accepted(record,1)=slice.fcn_lik(z_old_R);


    else
        Sigmadraw=SigmaHat;

    end

    if draw_B==1
        %% Draw B
        %cholSigmadraw = hh(Sigmadraw)';
        %Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*nvar,1) + reshape(mmuTilde,nvar*m,1);
        %Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);

        cholLSigmadraw = hh(Sigmadraw)';
        %cholLvarB = kron(cholLSigmadraw,cholOomegaTilde);
        slice.scale_z    = 1;
        slice.mean       = vec(mmuTilde);
        %slice.chol_cov_z = cholLvarB;%chol(Rvariance,'lower');
		slice.chol_cov_z1 = cholLSigmadraw;%chol(Rvariance,'lower');       
        slice.chol_cov_z2 = cholOomegaTilde;%chol(Rvariance,'lower');
        slice.nobs1 = info.nvar;
        slice.nobs2 = info.m;																   

        slice.fcn_lik    = @(z_prop) loglike_B(z_prop,function_restrictions_i,fo_inv,fo_str2irfs,Sigmadraw,Qdraw,info,mmuTilde,OomegaTilde);

        slice.nobs       = info.nvar*info.m;

        % --- actual slice sampling

        lik_old = slice.fcn_lik(z_old_B);

        %[z_old_B, ~, n_try] = slice_sampling_v02(slice, z_old_B, lik_old);
		[z_old_B, ~, n_try] = slice_sampling_v02_kron(slice, z_old_B, lik_old);																	

        Bdraw = reshape(z_old_B,info.m,info.nvar);

    else
        Bdraw         = BHat ;
    end


    x            = [vec(Bdraw); vec(Sigmadraw); vec(Qdraw)];
    structpara   = f_h_inv(x,info);
    irfpara      = fo_str2irfs(structpara);


    % if counter_Q==1000
    %        counter_Q=0;
    % end



    %% check if sign restrictions hold

    L0 =reshape(irfpara(1:nvar*nvar,:),nvar,nvar);

    BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];




    count=count+1;


        if rem(record, save_every) == 0
        % --- store
    % store orthogonal reduced-form draws
       ind_save = ind_save + 1;
    Bdraws{ind_save,1}     = Bdraw;
    Sigmadraws{ind_save,1} = Sigmadraw;
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
L    = zeros(horizon+1,nvar,nvar,ind_save);
cumL    = zeros(horizon+1,nvar,nvar,ind_save);


for s=1:ind_save

 

        Bdraw =     Bdraws{s,1} ;
        Sigmadraw = Sigmadraws{s,1} ;
        Qdraw=Qdraws{s,1};
 
        BSigmaQ = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];



    structpara = f_h_inv(BSigmaQ,info);


    LIRF = IRF_horizons(structpara, nvar, nlag, m, nex, 0:horizon);


 


    for h=0:horizon
        L(h+1,:,:,s) =  LIRF(1+h*nvar:(h+1)*nvar,:);
      
        for i=1:nvar
            cumL(h+1,1:4,i,s)   = sum(L(1:h+1,1:4,i,s),1);
        end

    end



end

stack = cat(3, Sigmadraws{:});        % 10×3×3 array
SigmaMean = mean(stack, 3);  

stack = cat(3, Bdraws{:}); 
BMean = mean(stack, 3);  


% cd results
% 
%    if ~exist('matfiles', 'dir')
%        mkdir('matfiles')
%    end
% 
% 
%     savefile= ['matfiles/rgibbs_',label_R,'prior_only_',num2str(prior_only),'prior_',prior_type,'fixed_rf_',num2str(fixed_rf),'_ndraws',num2str(M0),'n_shocks',num2str(info.nshocks),'.mat'];
% 
% save(savefile,'L','cumL','horizon','telapsed','-v7.3');
% cd ..

%% Efficiency factor
h=1;
% i=3;
% j=1;
ineff = zeros(size(L(:,:,:,1)));
for h = 1:1
    for i = 1:size(L,2)
        for j = 1:size(L,3)
            [a,b] = autocorr(squeeze(L(h,i,j,500:1:end)), 'NumLags',100);
            ineff(h,i,j) = 1 + 2*sum(a(2:end));
        end
    end
end

disp('---')
median(vec(squeeze(ineff(1,:,1:info.nshocks))))
quantile(vec(squeeze(ineff(1,:,1:info.nshocks))), 0.95)
quantile(vec(squeeze(ineff(1,:,1:info.nshocks))), 0.05)


plot(vec(squeeze(ineff(1,:,1:info.nshocks))))