function y = negative_loglike_obj_mp(x,info,Xnd,Ynd,MIN,MAX)


llambda0 = MIN.llambda+(MAX.llambda-MIN.llambda)./(1+exp(-x(1,1)));
pphi0    = MIN.pphi+(MAX.pphi-MIN.pphi)./(1+exp(-x(2:info.nvar+1,1)));




% set position where variables are in first differences


% inverse Wishart prior parameters
nnuBar              = info.nvar + 2;
PphiBar             = diag(pphi0)*(nnuBar-info.nvar-1);
oomegaBar=zeros(info.m,1);
oomegaBar(1:info.nex,1)=info.Vc;
for ell=1:info.nlag
    oomegaBar(info.nex+1+(ell-1)*info.nvar:info.nex+ell*info.nvar,1) =  (llambda0^2)./((ell^2)*pphi0);
end

OomegaBar=diag(oomegaBar);
OomegaBarInverse=diag(1./oomegaBar);
PpsiBar = zeros(info.m,info.nvar);
diagPpsiBar=ones(info.nvar,1);
diagPpsiBar(info.pos)=0;                % Set to zero the prior mean on the first own lag for variables selected in the vector pos
PpsiBar(info.nex+1:info.nvar+info.nex,:)=diag(diagPpsiBar);



% compute the marginal likelihood
Y = Ynd;
X = Xnd;



T = size(Y,1);

nnuTilde            = T  + nnuBar;
OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(info.m);
OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
PpsiTilde           = (X'*X  + OomegaBarInverse)\(X'*Y + OomegaBarInverse*PpsiBar);


PphiTilde           = Y'*Y + PphiBar + PpsiBar'*OomegaBarInverse*PpsiBar - PpsiTilde'*OomegaTildeInverse*PpsiTilde;
PphiTilde           = (PphiTilde'+PphiTilde)*0.5;

%tmp1 = -info.nvar*T*0.5*log(2*pi) + cIWln(nnuBar,PphiBar)-cIWln(nnuTilde,PphiTilde);

%tmp2 = -info.nvar*0.5*LogAbsDet(OomegaBar) + info.nvar*0.5*LogAbsDet(OomegaTilde);


% numerical stability tricks
eeps = Y-X*PpsiTilde;
DPhi=chol(PphiBar\eye(info.nvar),'lower');
%tmp1_stable = -T*0.5*LogAbsDet(PphiBar) - (T+nnuBar)*0.5*LogAbsDet(eye(info.nvar)+DPhi'*(eeps'*eeps  + (PpsiTilde-PpsiBar)'*OomegaBarInverse*(PpsiTilde-PpsiBar))*DPhi);
bbb=DPhi'*(eeps'*eeps  + (PpsiTilde-PpsiBar)'*OomegaBarInverse*(PpsiTilde-PpsiBar))*DPhi;
eigbbb=real(eig(bbb));
eigbbb(eigbbb<1e-12)=0;
eigbbb=eigbbb+1;
tmp1_stable = -T*0.5*LogAbsDet(PphiBar)- (T+nnuBar)*0.5*sum(log(eigbbb));

D=chol(OomegaBar,'lower');
%tmp2_stable = -info.nvar*0.5*LogAbsDet(D'*(X'*X)*D + eye(info.m));
aaa=D'*(X'*X)*D;
eigaaa=real(eig(aaa));
eigaaa(eigaaa<1e-12)=0;
eigaaa=eigaaa+1;
tmp2_stable =  -info.nvar*0.5*sum(log(eigaaa));

tmpg = 0;
for i=1:info.nvar
    tmpg =  tmpg  + gammaln((T+nnuBar+1-i)/2)-gammaln((nnuBar+1-i)/2);
end

yextended_stable = -info.nvar*T*0.5*log(pi) + tmpg + tmp1_stable + tmp2_stable;


y = (yextended_stable);

% hyperparameters
y=y+logGammapdf(llambda0,info.priorcoef.llambda.k,info.priorcoef.llambda.theta);
y=y+sum(logIG2pdf(pphi0/(nnuBar-info.nvar-1),info.priorcoef.alpha.PHI,info.priorcoef.beta.PHI));
y=y*-1;

function r=logGammapdf(x,k,theta);
r=(k-1)*log(x)-x/theta-k*log(theta)-gammaln(k);

function r=logIG2pdf(x,alpha,beta);
r=alpha*log(beta)-(alpha+1)*log(x)-beta./x-gammaln(alpha);