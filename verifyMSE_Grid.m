%This file is created by Xu Xiaoli on 22/06/2022
%It verifies the AMSE for channel prediction with PPP sample locations. 

clc;
clear;
close all;

%Specify the channel parameters 
alpha=8; %shadowing power
beta=30; %shadowing correlation distance
n_PL=2.2; %path loss exponent
sigma=sqrt(2); %multipath variance
KdB=-80; %reference received power
theta=[KdB; n_PL];
%D=500;
d_vec=[5 10 25];
k_vec=1:10;
numSamples_lb=1000; % the minimum number of samples in the CGM

numTarget=500;
AMSE_simu=zeros(length(d_vec),length(k_vec)+1); %also include k=0
AMSE_simu_formula=zeros(length(d_vec),length(k_vec)+1);
iter=10;
for p=1:iter
    for j=1:length(d_vec)
        d=d_vec(j); %the minimum separateion
        D=sqrt(numSamples_lb)*d; % determine the region size to ensure there are at least 10 times the samples 
        MeasureLoc=getMeasureLoc_Grid(D,d);
        numSamples=length(MeasureLoc(1,:));%the total number of sample points
        Totalnum=numSamples+numTarget;
        TargetLocation=[(rand(1,numTarget)-1/2)*D;(rand(1,numTarget)-1/2)*D]; 
        AllPointLoc=[MeasureLoc,TargetLocation];

        %Generate the values for those locations
        distance=sqrt(AllPointLoc(1,:).^2+AllPointLoc(2,:).^2); %BS at origin
        Hq_all=[ones(Totalnum,1),-10*log10(distance')];
        Hq_sample=[ones(numSamples,1),-10*log10(distance(1:numSamples)')];
        distMatrixAll=sqrt((repmat(AllPointLoc(1,:),Totalnum,1)-repmat(AllPointLoc(1,:)',1,Totalnum)).^2....
            +(repmat(AllPointLoc(2,:),Totalnum,1)-repmat(AllPointLoc(2,:)',1,Totalnum)).^2); %the distance matrix for both samples and targets
        RQ_all=alpha*exp(-distMatrixAll/beta);
        YQ_all=Hq_all*theta+(mvnrnd(zeros(1,Totalnum),RQ_all,1))'...
            +(mvnrnd(zeros(1,Totalnum),sigma^2*eye(Totalnum),1))';
        YQ_sample=YQ_all(1:numSamples);
        YQ_target=YQ_all(numSamples+1:end);

        %===Simulate the error the each target
        EstiPower=zeros(numTarget,length(k_vec));
        MSE_formula=zeros(numTarget,length(k_vec));
        for i=1:numTarget
            q_Loc=TargetLocation(:,i);
            hq=Hq_all(numSamples+i,:);
            %Estimated power for k=1
            EstiPower(i,1)=hq*theta;
            MSE_formula(i,1)=alpha+sigma^2;

            %Estimated power for larger k
            for u=1:length(k_vec)
                k=k_vec(u);
                [distk,distMatrix,~,idx]=findkNearest(q_Loc,MeasureLoc,k);
                Hq=Hq_sample(idx,:);
                PhiQ=alpha*exp(-distk'/beta);
                RQ=alpha*exp(-distMatrix/beta)+sigma^2*eye(k);
                tmp=PhiQ'/RQ;
                EstiPower(i,u+1)=EstiPower(i,1)+tmp*(YQ_sample(idx)-Hq*theta);
                MSE_formula(i,u+1)= MSE_formula(i,1)-tmp*PhiQ;
            end
        end
        AMSE_simu(j,:)=AMSE_simu(j,:)+sum((EstiPower-repmat(YQ_target,1,length(k_vec)+1)).^2)/numTarget;
        AMSE_simu_formula(j,:)=AMSE_simu_formula(j,:)+sum(MSE_formula)/numTarget;
    end
end
AMSE_simu=AMSE_simu/iter;
AMSE_simu_formula=AMSE_simu_formula/iter;
%Get the analysis results



figure;
plot([0,k_vec],AMSE_simu_formula(1,:),'rs-');
hold on;
plot([0,k_vec],MSEGridana(d_vec(1),[0,k_vec],alpha,beta,sigma),'r--');

plot([0,k_vec],AMSE_simu_formula(2,:),'bo-');
plot([0,k_vec],MSEGridana(d_vec(2),[0,k_vec],alpha,beta,sigma),'b--');
plot([0,k_vec],AMSE_simu_formula(3,:),'mv-');
plot([0,k_vec],MSEGridana(d_vec(3),[0,k_vec],alpha,beta,sigma),'m--');
hold off;
xlabel('$k$');
ylabel('AMSE');
legend('$d=5$ (simu)','$d=5$  (ana)','$d=10$ (simu)','$d=10$ (ana)','$d=25$ (simu)','$d=25$ (ana)','interpreter','latex')
grid on;
%Specify the number of test points
save GridMSE.mat;



