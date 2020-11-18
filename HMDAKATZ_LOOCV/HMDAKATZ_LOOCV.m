function HMDAKATZ_LOOCV(gamadd,gamall,beta,k)
%predict drug-microbe associations based on KATZ in the term of 5-fold cross validation
%A: Binary relations between drug and microbe, 1st column:drug, 2nd column:microbe

load knowndrugmicrobeinteraction.mat;
A=dd;
% nd:the number of drugs
% nm:the number of microbe
% pp:the number of known microbe-drug associations
nd=max(A(:,1)); 
nm=max(A(:,2));
[pp,qq]=size(A);
%interaction: adjacency matrix for the microbe-drug association network
%interaction(i,j)=1 means microbe j is related to drug i
for i=1:pp
    interaction(A(i,1),A(i,2))=1;
end

save interaction interaction;

load drugsmiles.mat;
     
%implement leave-one-out cross validation
for cv=1:pp 
    % obtain training sample
    load interaction.mat;
    interaction(A(cv,1),A(cv,2))=0;
   
    %calculate gamad for Gaussian kernel calculation
    for i=1:nd
        sd(i)=norm(interaction(i,:))^2;
    end
    gamad=nd/sum(sd')*gamadd;
    
    %calculate gamal for Gaussian kernel calculation
    for i=1:nm
        sl(i)=norm(interaction(:,i))^2;
    end
    gamal=nm/sum(sl')*gamall;
    
    %calculate Gaussian kernel for the similarity between drug: kd
    for i=1:nd
        for j=1:nd
            kd(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);
        end
    end
    
    %calculate Gaussian kernel for the similarity between microbe: km
    for i=1:nm
        for j=1:nm
            km(i,j)=exp(-gamal*(norm(interaction(:,i)-interaction(:,j)))^2);
        end
    end 

    kd = (kd +drugsmiles)/2;
        
    if k==2
        F=beta*interaction'+(beta^2)*(km*interaction'+interaction'*kd);
    else
        if k==3
            F=beta*interaction'+(beta^2)*(km*interaction'+interaction'*kd)+(beta^3)*(interaction'*interaction*interaction'+km*km*interaction'+km*interaction'*kd+interaction'*kd*kd);
        else
            if k==4
                 F=beta*interaction'+(beta^2)*(km*interaction'+interaction'*kd)+(beta^3)*(interaction'*interaction*interaction'+km*km*interaction'+km*interaction'*kd+interaction'*kd*kd)+(beta^4)*(km*km*km*interaction'+interaction'*interaction*km*interaction'+km*interaction'*interaction*interaction'+interaction'*kd*interaction*interaction')+(beta^4)*(interaction'*interaction*interaction'*kd+km*km*interaction'*kd+km*interaction'*kd*kd+interaction'*kd*kd*kd);
            end
        end
    end
    

% obtain the score of tested microbe-drug associations
finalscore=F(A(cv,2),A(cv,1));
% make the score of seed microbe-drug associations as zero
for i=1:nd
    for j=1:nm
        if interaction(i,j)==1
           F(j,i)=-10000;
        end
    end
end

% obtain the position of tested microbe-drug association as variable globalposition(1,cv),
[ll1,mm1]=size(find(F>=finalscore));
[ll2,mm2]=size(find(F>finalscore));
globalposition(1,cv)=ll2+1+(ll1-ll2-1)/2;

end
save('globalposition.mat','globalposition');   
end                