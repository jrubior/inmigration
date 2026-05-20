function loglik_prop = loglike_ttheta(z_prop,function_restrictions,fo_inv,fo_str2irfs,ttheta_old,sigma2draw,Qdraw,info,posterior_redu)
% likelihood eval y given z_prop
% z_prop = [ddelta, ggamma]


% ---


break0=0;
nvar=size(sigma2draw,1);
for i=1:nvar
            break1=break0+posterior_redu.ki{i};
            ttheta_old{i}=z_prop(break0+1:break1,1);
            break0=break1;
end

[Bdraw,Sigmadraw] = tthetassigma2TOBSigma(ttheta_old,sigma2draw,info);

S_prop = function_restrictions(Bdraw, Sigmadraw, Qdraw,fo_inv,fo_str2irfs,info);




% p(B|Sigma)
%varBgivenSigma = kron(Sigmadraw,OomegaTilde);
%varBgivenSigma = (varBgivenSigma+varBgivenSigma')/2;
 
%logpdfBgivenSigma = logmvnpdf_mc(z_prop,vec(mmuTilde),varBgivenSigma);%log(mvnpdf(vec(Bdraw),vec(mmuTilde),varBgivenSigma));


loglik_prop = -Inf;
if S_prop==1
loglik_prop = log(1);%logpdfBgivenSigma ;
end
