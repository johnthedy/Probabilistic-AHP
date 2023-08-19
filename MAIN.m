clc;clear
%% Input
%Consistency Index
RIlist(:,1)=1:1:10;
RIlist(:,2)=[0 0 0.58 0.9 1.12 1.24 1.32 1.41 1.45 1.49];

%Accepted Consistency Ratio
CRlimit=0.1;

%Choose problem, 3 problem provided in problem.m
nprob=1;
[nAHP,AHPraw]=problem(nprob);

%% Define Beta distribution
AHPdistparameter=AHPraw;
ULbound=[];
LLbound=[];       
for i=1:nAHP-1
    for j=i+1:nAHP
        dummy=AHPraw{i,j};
        if length(dummy)==1
            AHPdistparameter{i,j}=AHPraw{i,j};
        else
            for k=1:length(dummy)
                if dummy(k) >= 1
                    W(k) = 10-dummy(k);
                elseif dummy(k) < 1
                    W(k) = 1/dummy(k) + 8;
                end
            end
            mu=mean(W);
            varr=var(W);
            UL=max(W);
            LL=min(W);
            alpha=(mu-LL)/(UL-LL)*(((mu-LL)/(UL-LL))*(1-((mu-LL)/(UL-LL)))/(varr/((UL-LL)^2))-1);
            beta=(1-(mu-LL)/(UL-LL))*(((mu-LL)/(UL-LL))*(1-((mu-LL)/(UL-LL)))/(varr/((UL-LL)^2))-1);
            AHPdistparameter{i,j}=[alpha beta 1];
            LLbound=[LLbound,LL];
            ULbound=[ULbound,UL];
            clear W
        end
    end
end

clear dummy mu varr

%% Perform SOS optimization
ite=300;
ecosize=100;

n=length(LLbound);
ub=ULbound-0.1;
lb=LLbound+0.1;

eco=rand(ecosize,n).*(ub-lb)+lb;
fitness=zeros(ecosize,1);
for i=1:ecosize
    fitness(i,:)=fobj(eco(i,:),nAHP,LLbound,ULbound,AHPdistparameter,RIlist,CRlimit); 
end

for h=1:ite
    for i=1:ecosize
        % Update the best Organism
        [bestFitness,idx]=min(fitness); bestOrganism=eco(idx,:);
        
        %Mutualism Phase
        j=i;
        while i==j
            seed=randperm(ecosize); 
            j=seed(1);                  
        end
        % Determine Mutual Vector & Beneficial Factor
        mutualVector=mean([eco(i,:);eco(j,:)]);
        BF1=round(1+rand); BF2=round(1+rand);
        % Calculate new solution after Mutualism Phase
        ecoNew1=eco(i,:)+rand(1,n).*(bestOrganism-BF1.*mutualVector); 
        ecoNew2=eco(j,:)+rand(1,n).*(bestOrganism-BF2.*mutualVector);
        ecoNew1=bound(ecoNew1,ub,lb); 
        ecoNew2=bound(ecoNew2,ub,lb);
        % Evaluate the fitness of the new solution
        fitnessNew1=fobj(ecoNew1,nAHP,LLbound,ULbound,AHPdistparameter,RIlist,CRlimit);
        fitnessNew2=fobj(ecoNew2,nAHP,LLbound,ULbound,AHPdistparameter,RIlist,CRlimit);
        % Accept the new solution if the fitness is better
        if fitnessNew1<fitness(i)
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
        end
        if fitnessNew2<fitness(j)
           fitness(j)=fitnessNew2;
           eco(j,:)=ecoNew2;
        end
        % End of Mutualism Phase 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        % Commensialism Phase
        j=i;
        while i==j
            seed=randperm(ecosize); 
            j=seed(1);                  
        end
        % Calculate new solution after Commensalism Phase    
        ecoNew1=eco(i,:)+(rand(1,n)*2-1).*(bestOrganism-eco(j,:));
        ecoNew1=bound(ecoNew1,ub,lb);
        % Evaluate the fitness of the new solution
        fitnessNew1=fobj(ecoNew1,nAHP,LLbound,ULbound,AHPdistparameter,RIlist,CRlimit);
        % Accept the new solution if the fitness is better
        if fitnessNew1<fitness(i)
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
        end
        % End of Commensalism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Parasitism Phase
        j=i;
        while i==j
            seed=randperm(ecosize);
            j=seed(1);
        end
        % Determine Parasite Vector & Calculate the fitness
        parasiteVector=eco(i,:);
        seed=randperm(n);           
        pick=seed(1:ceil(rand*n));  % select random dimension
        parasiteVector(:,pick)=rand(1,length(pick)).*(ub(pick)-lb(pick))+lb(pick);
        fitnessParasite=fobj(parasiteVector,nAHP,LLbound,ULbound,AHPdistparameter,RIlist,CRlimit);

        % Kill organism j and replace it with the parasite 
        % if the fitness is lower than the parasite
        if fitnessParasite < fitness(j)
            fitness(j)=fitnessParasite;
            eco(j,:)=parasiteVector;
        end
        % End of Parasitism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
    end
    bestFitnesshistory(h)=bestFitness;
end
[bestFitness,idx]=min(fitness); bestOrganism=eco(idx,:);

AHPvariable=bestOrganism;
AHP=eye(nAHP);index=1;
for i=1:nAHP-1
    for j=i+1:nAHP
        if length(AHPdistparameter{i,j})==3
            dummy=AHPdistparameter{i,j};
            if AHPvariable(index)<=9
                fitness=fitness+betapdf((AHPvariable(index)-LLbound(index))/(ULbound(index)-LLbound(index)),dummy(1),dummy(2));
                betapdf((AHPvariable(index)-LLbound(index))/(ULbound(index)-LLbound(index)),dummy(1),dummy(2));
                AHP(i,j)=10-AHPvariable(index);
                AHP(j,i)=1/AHP(i,j);
            else
                fitness=fitness+betapdf((AHPvariable(index)-LLbound(index))/(ULbound(index)-LLbound(index)),dummy(1),dummy(2));
                betapdf((AHPvariable(index)-LLbound(index))/(ULbound(index)-LLbound(index)),dummy(1),dummy(2));
                AHP(i,j)=1/(AHPvariable(index)-8);
                AHP(j,i)=1/AHP(i,j);
            end
            index=index+1;
        end
    end
end

[eigenvector,eigenvalue]=eig(AHP);
[maxeigenvalue,maxorder]=max(max(abs(diag(eigenvalue))));
CI=(maxeigenvalue-nAHP)/(nAHP-1);
RI=RIlist(RIlist(:,1)==nAHP,2);
CR=CI/RI;
ratio=abs(eigenvector(:,maxorder))./sum(abs(eigenvector(:,maxorder)));

%% Ouput
disp(ratio)
disp(CR)