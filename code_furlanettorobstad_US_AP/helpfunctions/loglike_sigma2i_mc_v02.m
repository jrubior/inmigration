function loglik_prop = loglike_sigma2i_mc_v02(z_prop,function_restrictions,fo_inv,fo_str2irfs,ttheta_old,sigma2draw,Qdraw,info,nnuTilde,ii)
% likelihood eval y given z_prop
% z_prop = [ddelta, ggamma]


%[Bdraw,Sigmadraw] = tthetassigma2TOBSigma(ttheta_old,sigma2draw,info);

sigma2draw(ii,1)=1/(z_prop'*z_prop);


[Bdraw,Sigmadraw] = tthetassigma2TOBSigma(ttheta_old,sigma2draw,info);

if min(eig(Sigmadraw))>1e-8

S_prop = function_restrictions(Bdraw, Sigmadraw, Qdraw,fo_inv,fo_str2irfs,info);





loglik_prop = -Inf;
if S_prop==1 
    
loglik_prop = 0;
end

else
    loglik_prop = -Inf;

end

