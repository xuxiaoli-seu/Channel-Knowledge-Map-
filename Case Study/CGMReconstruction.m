%This file is created by Xu Xiaoli on 19/10/2023
%It reconstruct the CGM use Model-based channel prediciton 

clc;
clear;
close all;

load MeasureData2_4GHz_d20.mat;

numSamples=length(MeasureLoc(1,:));
distance=sqrt(MeasureLoc(1,:).^2+MeasureLoc(2,:).^2);
Hq_sample=[ones(numSamples,1),-10*log10(distance')];

%====Path loss based on the whole area======
LoSIdx=find(LoS==1); %The index for LoS points
NLoSIdx=setdiff(1:numSamples,LoSIdx); %The index for NLoS points

Hq_LoS=Hq_sample(LoSIdx,:);
YQ_LoS=YQ(LoSIdx);
MeasureLoc_LoS=MeasureLoc(:,LoSIdx);

tmp=Hq_LoS/(Hq_LoS'*Hq_LoS);
theta_LoS=YQ_LoS*tmp;
YQ_LoS_sha=YQ_LoS'-(Hq_LoS*theta_LoS');
%[alpha_LoS, beta_LoS,sigmasq_LoS]=EstShadowPara_LS(YQ_LoS_sha,MeasureLoc_LoS); %estimate the shadowing power
%no shadowing results alpha_LoS=0, beta_LoS=0

Hq_NLoS=Hq_sample(NLoSIdx,:);
YQ_NLoS=YQ(NLoSIdx);
MeasureLoc_NLoS=MeasureLoc(:,NLoSIdx);

tmp=Hq_NLoS/(Hq_NLoS'*Hq_NLoS);
theta_NLoS=YQ_NLoS*tmp;
YQ_NLoS_sha=YQ_NLoS'-(Hq_NLoS*theta_NLoS');
[alpha_NLoS, beta_NLoS,sigmasq_NLoS]=EstShadowPara_LS(YQ_NLoS_sha,MeasureLoc_NLoS); %estimate the shadowing power

%theta_est =[ -15.9144    5.2013]; 
figure;
plot(MeasureLoc_LoS(1,:),MeasureLoc_LoS(2,:),'rs','MarkerFaceColor','r','MarkerSize',2);
hold on;
plot(MeasureLoc_NLoS(1,:),MeasureLoc_NLoS(2,:),'bo','MarkerFaceColor','b','MarkerSize',2);
plot(0,0,'rp','MarkerFaceColor','r');
xlim([-140,250]);
ylim([-165, 250]);

%=====Neighboring region======
NumRemain=length(YQ_remain); %the number of points that need to be reconstructed 
k=3; %the neighbors used for prediction
YQ_estimated=zeros(1,NumRemain);
for i=1:NumRemain
    Loc=MeasureLoc_remain(:,i);
    hq=[1,-10*log10(norm(Loc))];
    if LoS_remain(i)==1 %this is LoS point
        YQ_estimated(i)=hq*theta_LoS';
    else
        [distk,distMatrix,~,idx]=findkNearest(Loc,MeasureLoc_NLoS,k);
        Hq=Hq_NLoS(idx,:);
        PhiQ=alpha_NLoS*exp(-distk'/beta_NLoS);
        RQ=alpha_NLoS*exp(-distMatrix/beta_NLoS);
        tmp=PhiQ'/RQ;
        YQ_estimated(i)=hq*theta_NLoS'+tmp*(YQ_NLoS(idx)'-Hq*theta_NLoS');
    end
end

AMSE_simu=sum((YQ_estimated-YQ_remain).^2)/NumRemain

%Substitute the estiamted value back and view the results
YQ_remain_all(Indx_keep)=YQ_estimated;

remainMeasurement=reshape(YQ_remain_all,row,col);
ChannelGainMatrix(remainPointIdx,remainPointIdx)=remainMeasurement';
XLoc=XRange-TxLoc(1);
YLoc=YRange-TxLoc(2);
[X,Y]=meshgrid(XLoc,YLoc);
%==========View the LoS
s=surf(X,Y,ChannelGainMatrix');
s.EdgeColor = 'none';
colorbar;
hold on;

view(2)
xlim([-150,280]);
ylim([-180,250]);
plot(0,0,'rp','MarkerFaceColor','r');
