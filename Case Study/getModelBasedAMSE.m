%This file is created by Xu Xiaoli on 2/11/2023
%It returns the AMSE give a set of sample Loc and Test Loc

function AMSE=getModelBasedAMSE(SampleLoc,YQ_sample,TestLoc,YQ_test)

numSamples=length(YQ_sample);
distance=sqrt(SampleLoc(1,:).^2+SampleLoc(2,:).^2);
Hq_sample=[ones(numSamples,1),-10*log10(distance')];
tmp=Hq_sample/(Hq_sample'*Hq_sample);
theta=YQ_sample*tmp;
if theta(2)<0
    %negative path loss
    theta=[sum(YQ_sample)/numSamples, 0];
end
theta
YQ_sha=YQ_sample'-(Hq_sample*theta');
[alpha, beta,sigmasq]=EstShadowPara_LS(YQ_sha,SampleLoc) 

NumRemain=length(YQ_test); %the number of points that need to be reconstructed 
k=3; %the neighbors used for prediction
YQ_estimated=zeros(1,NumRemain);
for i=1:NumRemain
    Loc=TestLoc(:,i);
    hq=[1,-10*log10(norm(Loc))];
    if alpha<0.01 %this is LoS point
        YQ_estimated(i)=hq*theta';
    else
        [distk,distMatrix,~,idx]=findkNearest(Loc,SampleLoc,k);
        Hq=Hq_sample(idx,:);
        PhiQ=alpha*exp(-distk'/beta);
        RQ=alpha*exp(-distMatrix/beta);
        tmp=PhiQ'/RQ;
        YQ_estimated(i)=hq*theta'+tmp*(YQ_sample(idx)'-Hq*theta');
    end
end

AMSE=sum((YQ_estimated-YQ_test).^2)/NumRemain;