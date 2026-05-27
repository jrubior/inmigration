function [loglik_prop,Bdraw,Sigmadraw,Adraw,Bstrdraw] = loglike_sigma2i_mc_BSig_v02(z_prop,function_restrictions,ttheta_old,sigma2draw,Qdraw,info,nnuTilde,ii, Bold, Sigold, Aold,Bstrold)
% likelihood eval y given z_prop
% z_prop = [ddelta, ggamma]


%[Bdraw,Sigmadraw] = tthetassigma2TOBSigma(ttheta_old,sigma2draw,info);

sigma2draw(ii,1)=1/(z_prop'*z_prop);


[Bdraw,Sigmadraw,Adraw,Bstrdraw] = tthetassigma2TOBSigma_BSig(ttheta_old,sigma2draw,info,ii,Bold,Sigold,Aold,Bstrold);

if min(eig(Sigmadraw))>1e-8
    
    S_prop = function_restrictions(Bdraw, Sigmadraw, Qdraw,info);
    
    loglik_prop = -Inf;
    if S_prop==1
        

        loglik_prop = -0.5*LogAbsDet(sigma2draw(ii,1)*posterior_redu.Kthetai{ii}\eye(posterior_redu.ki{ii}))-0.5*(ttheta_old{ii}-posterior_redu.thetai_tilde{ii})'*(posterior_redu.Kthetai{ii}/sigma2draw(ii,1))*(ttheta_old{ii}-posterior_redu.thetai_tilde{ii});


        %loglik_prop = 0;
    end
    
else
    loglik_prop = -Inf;
    
end

