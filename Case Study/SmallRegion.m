%This file is created by Xu Xiaoli on 2/11/2023
%It reconstruct the CGM use Model-based channel prediciton 
%It considers the reconstruction with a small area


clc;
clear;
close all;

load MeasureData2_4GHz_d20.mat;

%The samples wiithin the interested region
Indx=find((MeasureLoc(1,:)<60 ).*(MeasureLoc(1,:)>-60).*(MeasureLoc(2,:)>-90).*(MeasureLoc(2,:)<0));
SampleLoc=MeasureLoc(:,Indx);
YQ_sample=YQ(Indx);
%The test point within the interested region
Indx2=find((MeasureLoc_remain(1,:)<60 ).*(MeasureLoc_remain(1,:)>-60).*(MeasureLoc_remain(2,:)>-90).*(MeasureLoc_remain(2,:)<0));
TestLoc=MeasureLoc_remain(:,Indx2);
YQ_test=YQ_remain(Indx2);

%Case I: no further division

AMSE=getModelBasedAMSE(SampleLoc,YQ_sample,TestLoc,YQ_test)
figure;
plot(SampleLoc(1,:),SampleLoc(2,:),'bo','MarkerFaceColor','b','MarkerSize',2);
hold on;

plot(0,0,'rp','MarkerFaceColor','r');


%===============================
%Divid into half plane
Indx1=find((SampleLoc(1,:)<0 )); %Group1
Indx2=setdiff(1:length(YQ_sample),Indx1); %Group 2
SampleLoc1=SampleLoc(:,Indx1);
SampleLoc2=SampleLoc(:,Indx2);
YQ_sample1=YQ_sample(Indx1);
YQ_sample2=YQ_sample(Indx2);
%Divide the test region
Indx1=find(TestLoc(1,:)<0);
Indx2=setdiff(1:length(YQ_test),Indx1); %Group 2
TestLoc1=TestLoc(:,Indx1);
TestLoc2=TestLoc(:,Indx2);
YQ_test1=YQ_test(Indx1);
YQ_test2=YQ_test(Indx2);

AMSE1=getModelBasedAMSE(SampleLoc1,YQ_sample1,TestLoc1,YQ_test1)
AMSE2=getModelBasedAMSE(SampleLoc2,YQ_sample2,TestLoc2,YQ_test2)

%===============================
%Divid into four parts
Indx1=find((SampleLoc1(2,:)<-50 )); %Group1
Indx2=setdiff(1:length(YQ_sample1),Indx1); %Group 2
SampleLoc11=SampleLoc1(:,Indx1);
SampleLoc12=SampleLoc1(:,Indx2);
YQ_sample11=YQ_sample1(Indx1);
YQ_sample12=YQ_sample1(Indx2);

Indx1=find((SampleLoc2(2,:)<-50 )); %Group1
Indx2=setdiff(1:length(YQ_sample2),Indx1); %Group 2
SampleLoc21=SampleLoc2(:,Indx1);
SampleLoc22=SampleLoc2(:,Indx2);
YQ_sample21=YQ_sample2(Indx1);
YQ_sample22=YQ_sample2(Indx2);

Indx1=find(TestLoc1(2,:)<-50);
Indx2=setdiff(1:length(YQ_test1),Indx1); %Group 2
TestLoc11=TestLoc1(:,Indx1);
TestLoc12=TestLoc1(:,Indx2);
YQ_test11=YQ_test1(Indx1);
YQ_test12=YQ_test1(Indx2);

Indx1=find(TestLoc2(2,:)<-50);
Indx2=setdiff(1:length(YQ_test2),Indx1); %Group 2
TestLoc21=TestLoc2(:,Indx1);
TestLoc22=TestLoc2(:,Indx2);
YQ_test21=YQ_test2(Indx1);
YQ_test22=YQ_test2(Indx2);

AMSE11=getModelBasedAMSE(SampleLoc11,YQ_sample11,TestLoc11,YQ_test11)
AMSE12=getModelBasedAMSE(SampleLoc12,YQ_sample12,TestLoc12,YQ_test12)
AMSE21=getModelBasedAMSE(SampleLoc21,YQ_sample21,TestLoc21,YQ_test21)
AMSE22=getModelBasedAMSE(SampleLoc22,YQ_sample22,TestLoc22,YQ_test22)
avgMSE=(AMSE11+AMSE12+AMSE21+AMSE22)/4

figure;
plot(SampleLoc(1,:),SampleLoc(2,:),'bs','MarkerFaceColor','b','MarkerSize',4);
hold on;
plot(TestLoc11(1,:),TestLoc11(2,:),'bo','MarkerFaceColor','b','MarkerSize',2);
plot(TestLoc12(1,:),TestLoc12(2,:),'ro','MarkerFaceColor','r','MarkerSize',2);
plot(TestLoc21(1,:),TestLoc21(2,:),'co','MarkerFaceColor','c','MarkerSize',2);
plot(TestLoc22(1,:),TestLoc22(2,:),'yo','MarkerFaceColor','y','MarkerSize',2);


plot(0,0,'rp','MarkerFaceColor','r');
xlim([-60,60]);
ylim([-90,0]);

figure;
x=["1 Region" "2 Regions" "4 Regions"];
X = categorical({'1 Region','2 Regions','4 Regions'});
X = reordercats(X,{'1 Region','2 Regions','4 Regions'});
%x=[1 2 3];
y=[0,AMSE,0,0;0, AMSE1,AMSE1,0; AMSE11,AMSE12,AMSE21,AMSE22];
b=bar(X,y);
hold on;
p1=plot(sum(y,2)./[1;2;4],'bo-');
p2=plot([109.9952,109.9952,109.9952],'k--');
legend([p1,p2],'Average AMSE','Area AMSE')
grid on;
return;

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
figure;
s=surf(X,Y,ChannelGainMatrix');
s.EdgeColor = 'none';
colorbar;
hold on;

view(2)
xlim([-150,280]);
ylim([-180,250]);
plot(0,0,'rp','MarkerFaceColor','r');
