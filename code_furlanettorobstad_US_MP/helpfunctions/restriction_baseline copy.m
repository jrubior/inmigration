function [St_old] = restriction_baseline(Bdraw, Sigmadraw, Qdraw,info)
    
% checking sign restrictions

% checking sign restrictions
nvar=info.nvar;


hSigmadraw = chol(Sigmadraw);

L0 = hSigmadraw'*Qdraw;



SRvec = [];
count = 1;


for i=1:info.nshocks

for j=1:rank(info.Ss{i})

SRvec(count,1) = info.Ss{i}(j,:)*L0*info.e(:,i);
count = count + 1;
end

end



St_old=all((SRvec>0)==1);

end








