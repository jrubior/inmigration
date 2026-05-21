function [nnuBar,PphiBar,PpsiBar,OomegaBarInverse,Y,X]  = SetupMinnesota(Ynd,Xnd,setup)  


% we follow the simplest prior in Giannone, Lenza, Primiceri (2015)
% minnesota prior only, that is no sum-of-coefficients and no dummy-initial observations

setup.pos = [];   % stationary variables (empty when all variables are assumed non-stationary)
setup.Vc  = 10^7; % variance of prior over the constant term c\sim N(0,Vc)

% Giannone, Lenza, Primiceri, (GLP,2015) propose to either incorporate
% hyperparameters as part of sampler or use MDD as likelihood for
% the hyper paraemeters and maximize posterior to obtain optimized
% hyper para that are treated as values. We adopt the latter.


% hyperparameters:
% lambda: governs overall tightness of prior for entries in B
% \Phi:   scale of covariance matrix, assumed diagonal: \phi =
% diag(\Phi), \Sigma ~ IW(\nu,\Phi)
% attention! GLP,2015 define the prior over pphi0/(nnuBar-info.nvar-1)


% Prior for hyperparameters:

% prior for \lambda is Gamma(k,theta) where k and theta are such that mode is 0.2 and std is 0.4.
mode.llambda = .2;
sd.llambda   = .4;
setup.priorcoef.llambda   = GammaCoef(mode.llambda,sd.llambda,0); % transform mode and std into Gamma(k,theta) coefficients


% prior for \Psi(i,i) is Gamma(alpha,beta) where \alpha = \beta = 0.02^2
scalePHI   = 0.02^2;
setup.priorcoef.alpha.PHI = scalePHI;
setup.priorcoef.beta.PHI  = scalePHI;


% Initial values for hyperparameters
llambda   = 0.5;   % initial guess for hyper-parameter \lambda


% initial guess for PHI = diag(phi): residual variance of AR(1) for each variable
SS=zeros(setup.nvar,1);
T = size(Ynd,1);
for i=1:setup.nvar
    yols=Ynd(2:end,i);
    xols=[ones(T-1,1),Ynd(1:end-1,i)];
    bbetaols = (xols'*xols)\(xols'*yols);
    eols = yols-xols*bbetaols;
    SS(i)=eols'*eols/(size(yols,1)-size(xols,2));
end
pphi      = SS;

% bounds for maximization
MIN.llambda = 0.0001;
MAX.llambda = 5;
MIN.pphi  = SS./100;
MAX.pphi  = SS.*100;


% transformation of initial values for optimization
inllambda = -log((MAX.llambda-llambda)./(llambda-MIN.llambda));
inppsi    = -log((MAX.pphi-pphi)./(pphi-MIN.pphi));

x0 = [inllambda;inppsi]; % vector to optimize


f0 = negative_loglike_obj_mp(x0,setup,Xnd,Ynd,MIN,MAX); % objective function: posterior \propto marginal data density is the likelihood times prior over hyper-parameters

crit                              = 1e-16;              % csminwel setting
nit                               = 1000;               % csminwel setting
H0                                = eye(size(x0,1))*10; % csminwel setting


if setup.reoptimizeMH==1

    [fh,xh,~,~,itct,~,~]              = csminwel('negative_loglike_obj_mp',x0,H0,[],crit,nit,setup,Xnd,Ynd,MIN,MAX);

    save('xh.mat','xh','fh');
else
    load('xh.mat','xh','fh');
end

llambda0_xh = MIN.llambda+(MAX.llambda-MIN.llambda)./(1+exp(-xh(1,1)));
pphi0_xh    = MIN.pphi+(MAX.pphi-MIN.pphi)./(1+exp(-xh(2:setup.nvar+1,1)));%

llambda0 =  llambda0_xh; % optimized lambda
pphi0    =  pphi0_xh;    % optimized phi


% we re-set the random seed after using csminwel to guarantee the
% posterior draws are identical across operating systems. Otherwise,
% csminwel affects the seed in an unpredictable manner across
% operating systems
rng(setup.seed,'twister');


% inverse Wishart prior parameters
nnuBar              = setup.nvar + 2;
PphiBar             = diag(pphi0);
oomegaBar=zeros(setup.m,1);
oomegaBar(1:setup.nex,1)=setup.Vc;
for ell=1:setup.nlag
    oomegaBar(setup.nex+1+(ell-1)*setup.nvar:setup.nex+ell*setup.nvar,1) =  (llambda0^2)*(nnuBar-setup.nvar-1)./((ell^2)*pphi0);
end
OomegaBar=diag(oomegaBar);
OomegaBarInverse=diag(1./oomegaBar);
PpsiBar = zeros(setup.m,setup.nvar);
diagPpsiBar=ones(setup.nvar,1);
diagPpsiBar(setup.pos)=0;   % Set to zero the prior mean on the first own lag for variables selected in the vector pos

PpsiBar(setup.nex+1:setup.nvar+setup.nex,:)=diag(diagPpsiBar);


Y = Ynd;
X = Xnd;


end