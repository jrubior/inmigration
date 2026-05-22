function [St_old] = restriction_baseline(Bdraw, Sigmadraw, Qdraw,setup)
    
% checking sign restrictions


hSigmadraw = chol(Sigmadraw);

L0 = hSigmadraw'*Qdraw;



SRvec = [];
count = 1;


for i=1:setup.nshocks

for j=1:rank(setup.Ss{i})

SRvec(count,1) = setup.Ss{i}(j,:)*L0*setup.e(:,i);
count = count + 1;
end

end



St_old=all((SRvec>0)==1);

end








