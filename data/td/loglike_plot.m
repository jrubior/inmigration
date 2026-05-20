function []=loglike_plot(res)

plot(res.r,res.loglik,'ro',res.r,res.loglik);
            title  ('Likelihood vs innovational/dynamic parameter');
            ylabel ('likelihood'); xlabel ('innovational/dynamic parameter');
            grid on;
            
            
