function HMDAKATZ_5CV(gamadd,gamall,beta,k)

%predict drug-microbe associations based on KATZ in the term of 5-fold cross validation
%A: Binary relations between drug and microbe, 1st column:drug, 2nd column:microbe

load knowndrugmicrobeinteraction.mat;
A=dd;
% nd:the number of drug
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

    %implement 5-fold cross validation
    x=randperm(pp)';
    T=1;

    N_fold = 5;
    for cv=1:N_fold
        load interaction interaction;
        if cv<N_fold
            B=A(x((cv-1)*floor(pp/N_fold)+1:floor(pp/N_fold)*cv),:);
            %obtain training sample
            for i=1:floor(pp/N_fold)
                interaction(B(i,1),B(i,2))=0;
            end
        else
            B=A(x((cv-1)*floor(pp/N_fold)+1:pp),:);
            %obtain training sample
            for i=1:pp-floor(pp/N_fold)*(N_fold-1)
                interaction(B(i,1),B(i,2))=0;
            end
        end

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
        
        kd=(kd+drugsmiles)/2;
        
        %default 
         Y=interaction;
            
        if k==2
            F=beta*Y'+(beta^2)*(km*Y'+Y'*kd);
        else
            if k==3
                F=beta*Y'+(beta^2)*(km*Y'+Y'*kd)+(beta^3)*(Y'*Y*Y'+km*km*Y'+km*Y'*kd+Y'*kd*kd);
            else
                if k==4
                    F=beta*Y'+(beta^2)*(km*Y'+Y'*kd)+(beta^3)*(Y'*Y*Y'+km*km*Y'+km*Y'*kd+Y'*kd*kd)+(beta^4)*(km*km*km*Y'+Y'*Y*km*Y'+km*Y'*Y*Y'+Y'*kd*Y*Y')+(beta^4)*(Y'*Y*Y'*kd+km*km*Y'*kd+km*Y'*kd*kd+Y'*kd*kd*kd);
                end
            end
        end
        F=F';
       
        [size1B,size2B]=size(B);
        % obtain the score of tested  microbe-drug associations 
        for i=1:size1B
            finalscore(i,1)=F(B(i,1),B(i,2));
        end
        %make the score of seed  microbe-drug associations as zero
        for i=1:nd
            for j=1:nm
                if interaction(i,j)==1
                   F(i,j)=-10000;
                end
            end
        end
        for qq=1:size1B
            %obtain the position of tested microbe-drug associations as variable position(1,cv), 
            [ll1,mm1]=size(find(F>=finalscore(qq)));
            [ll2,mm2]=size(find(F>finalscore(qq)));
            position(1,T)=ll2+1+(ll1-ll2-1)/2;
            T=T+1;
        end

    end
    save('position.mat','position');  

end
