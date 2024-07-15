%This file is created by Xu Xiaoli on 07/10/2023
%It return the simulated MSE based on the estimated parameters

function AMSE_simu=AMSEsimu_estiPara(MeasureLoc,TargetLocation, YQ_all,k,n_PL,KdB,alpha,beta)

% 
% global alpha; %shadowing power
% global beta; %shadowing correlation distance
% global n_PL; %path loss exponent
% global sigma; %multipath variance
% global KdB; %reference received power

theta=[KdB; n_PL];

numSamples=length(MeasureLoc);
numTarget=length(TargetLocation);
YQ_sample=YQ_all(1:numSamples);
YQ_target=YQ_all(numSamples+1:end);
AllPointLoc=[MeasureLoc,TargetLocation];
distance=sqrt(AllPointLoc(1,:).^2+AllPointLoc(2,:).^2); %BS at origin
Totalnum=numSamples+numTarget;
Hq_all=[ones(Totalnum,1),-10*log10(distance')];
Hq_sample=[ones(numSamples,1),-10*log10(distance(1:numSamples)')];

EstiPower=zeros(numTarget,1);
%MSE_ana=zeros(numTarget,1);
for i=1:numTarget
    q_Loc=TargetLocation(:,i);
    hq=Hq_all(numSamples+i,:);
    dist_vec=sqrt((MeasureLoc(1,:)-q_Loc(1)).^2+(MeasureLoc(2,:)-q_Loc(2)).^2);
    k=min(k,sum(dist_vec<10*beta));
    if k==0 
        EstiPower(i)=hq*theta; 
    else
       [distk,distMatrix,~,idx]=findkNearest(q_Loc,MeasureLoc,k);
        Hq=Hq_sample(idx,:);
        PhiQ=alpha*exp(-distk'/beta);
        %RQ=alpha*exp(-distMatrix/beta)+sigma^2*eye(k);
        RQ=alpha*exp(-distMatrix/beta);
        tmp=PhiQ'/RQ;
        EstiPower(i)=hq*theta+tmp*(YQ_sample(idx)-Hq*theta);
    end
    %MSE_ana(i)=alpha+sigma_sq-tmp*PhiQ;
    %
end

AMSE_simu=sum((EstiPower-YQ_target).^2)/numTarget;
%AMSE_ana=sum(MSE_ana)/numTarget;