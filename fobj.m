function fitness=fobj(AHPvariable,nAHP,LLbound,ULbound,AHPdistparameter,RIlist,CRlimit)

AHP=eye(nAHP);fitness=0;index=1;
for i=1:nAHP-1
    for j=i+1:nAHP
        if length(AHPdistparameter{i,j})==3
            dummy=AHPdistparameter{i,j};
            if AHPvariable(index)<=9
                fitness=fitness+betapdf((AHPvariable(index)-LLbound(index))/(ULbound(index)-LLbound(index)),dummy(1),dummy(2));
                AHP(i,j)=10-AHPvariable(index);
                AHP(j,i)=1/AHP(i,j);
            else
                fitness=fitness+betapdf((AHPvariable(index)-LLbound(index))/(ULbound(index)-LLbound(index)),dummy(1),dummy(2));
                AHP(i,j)=1/(AHPvariable(index)-8);
                AHP(j,i)=1/AHP(i,j);
            end
            index=index+1;
        end
    end
end
[~,eigenvalue]=eig(AHP);
[maxeigenvalue,~]=max(max(abs(diag(eigenvalue))));
CI=(maxeigenvalue-nAHP)/(nAHP-1);
RI=RIlist(RIlist(:,1)==nAHP,2);
CR=CI/RI;

if CR>CRlimit
    fitness=100;
else
    fitness=1/fitness;
end
end


