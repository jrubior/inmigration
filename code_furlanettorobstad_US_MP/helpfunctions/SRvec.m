function SRvec = SRvec(Bdraw,Sigmadraw,Qdraw,setup)


% checking sign restrictions
nvar=setup.nvar;


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









end