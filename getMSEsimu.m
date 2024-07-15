%This file is created by Xu Xiaoli on 22/06/2022
%It return the simulated MSE.

function [AMSE_simu, AMSE_ana]=getMSEsimu(MeasureLoc,numTarget,k)

global D;
global alpha; %shadowing power
global beta; %shadowing correlation distance
global n_PL; %path loss exponent
global sigma; %multipath variance
global KdB; %reference received power

theta=[KdB; n_PL];
numSamples=length(MeasureLoc(1,:));%the total number of sample points

Totalnum=numSamples+numTarget;
TargetLocation=[(rand(1,numTarget)-1/2)*D;(rand(1,numTarget)-1/2)*D]; 
AllPointLoc=[MeasureLoc,TargetLocation];

distance=sqrt(AllPointLoc(1,:).^2+AllPointLoc(2,:).^2); %BS at origin
Hq_all=[ones(Totalnum,1),-10*log10(distance')];
Hq_sample=[ones(numSamples,1),-10*log10(distance(1:numSamples)')];

distMatrixAll=sqrt((repmat(AllPointLoc(1,:),Totalnum,1)-repmat(AllPointLoc(1,:)',1,Totalnum)).^2....
    +(repmat(AllPointLoc(2,:),Totalnum,1)-repmat(AllPointLoc(2,:)',1,Totalnum)).^2); %the distance matrix for both samples and targets
%distanceMatrixSample=sqrt((repmat(MeasureLoc(1,:),numSamples,1)-repmat(MeasureLoc(1,:)',1,numSamples)).^2....
    %+(repmat(MeasureLoc(2,:),numSamples,1)-repmat(MeasureLoc(2,:)',1,numSamples)).^2); %the distance matrix only for the samples

RQ_all=alpha*exp(-distMatrixAll/beta);
%RQ_sample=alpha*exp(-distanceMatrixSample/beta);

YQ_all=Hq_all*theta+(mvnrnd(zeros(1,Totalnum),RQ_all,1))'...
    +(mvnrnd(zeros(1,Totalnum),sigma^2*eye(Totalnum),1))';
YQ_sample=YQ_all(1:numSamples);
YQ_target=YQ_all(numSamples+1:end);

EstiPower=zeros(numTarget,1);
MSE_ana=zeros(numTarget,1);
for i=1:numTarget
    q_Loc=TargetLocation(:,i);
    hq=Hq_all(numSamples+i,:);
    [distk,distMatrix,~,idx]=findkNearest(q_Loc,MeasureLoc,k);
    Hq=Hq_sample(idx,:);
    PhiQ=alpha*exp(-distk'/beta);
    %RQ=alpha*exp(-distMatrix/beta)+sigma^2*eye(k);
    RQ=alpha*exp(-distMatrix/beta);
    tmp=PhiQ'/RQ;
    EstiPower(i)=hq*theta+tmp*(YQ_sample(idx)-Hq*theta);
    
    MSE_ana(i)=alpha+sigma^2-tmp*PhiQ;
    %
end

AMSE_simu=sum((EstiPower-YQ_target).^2)/numTarget;
AMSE_ana=sum(MSE_ana)/numTarget;