function [loglik_prop,Bdraw,Sigmadraw,Adraw,Bstrdraw] = loglike_tthetai_BSig(z_prop,function_restrictions,fo_inv,fo_str2irfs,ttheta_old,sigma2draw,Qdraw,info,ki,ii,Bold,Sigold,Aold,Bstrold)
% likelihood eval y given z_prop
% z_prop = [ddelta, ggamma]


% ---



ttheta_old{ii}=z_prop;


[Bdraw,Sigmadraw,Adraw,Bstrdraw] = tthetassigma2TOBSigma_BSig(ttheta_old,sigma2draw,info,ii,Bold,Sigold,Aold,Bstrold);

S_prop = function_restrictions(Bdraw, Sigmadraw, Qdraw,fo_inv,fo_str2irfs,info);




% p(B|Sigma)
%varBgivenSigma = kron(Sigmadraw,OomegaTilde);
%varBgivenSigma = (varBgivenSigma+varBgivenSigma')/2;
 
%logpdfBgivenSigma = logmvnpdf_mc(z_prop,vec(mmuTilde),varBgivenSigma);%log(mvnpdf(vec(Bdraw),vec(mmuTilde),varBgivenSigma));


loglik_prop = -Inf;
if S_prop==1
loglik_prop = log(1);%logpdfBgivenSigma ;
end
