%This file is created by Xu Xiaoli on 14/07/2023
%It performs distributed Kriging for the CGM

clc;
clear;
close all;

load 'MeasureData2_4GHz_d5.mat';
distBS=sqrt(MeasureLoc(1,:).^2+MeasureLoc(2,:).^2);

%=================Select the number of test points=======================
numTestPoints=1000;
remainPoints=length(YQ_remain);
RandPermutation=randperm(remainPoints); %permutat the index
TestPointsIdx=RandPermutation(1:numTestPoints); %choose the set of test points
Loc_test=MeasureLoc_remain(:,TestPointsIdx); %those data not for trainning used for test
YQ_test=YQ_remain(TestPointsIdx);

%For each test point, use k neighboring measurements for estimating the
%gain

k_vec=5:10:200;

%for storing the estimated value with distributed Kriging
YQ_esti=zeros(length(k_vec),numTestPoint);
%for storing the average YQ
YQ_avg=zeros(length(k_vec),numTestPoint);

for i=1:numTestPoints
    q_Loc=Loc_test(:,i);
    %For complexity reduction, we constraint on a small region
    TrainRange=50;
    TrainIdx=(MeasureLoc(1,:)>q_Loc(1)-TrainRange)&(MeasureLoc(1,:)<q_Loc(1)+TrainRange)&(MeasureLoc(2,:)>q_Loc(2)-TrainRange)&(MeasureLoc(2,:)<q_Loc(2)+TrainRange);
    SampleLoc=MeasureLoc(:,TrainIdx);
    SampleYQ=YQ(TrainIdx); %The measurements at the sample
    
    for j=1:length(k_vec)
    %Get the k nearest neighbore
        k=k_vec(j);
        %====Find k nearest distance===============
        [distk,distMatrix,distBS,index]=findkNearest(q_Loc,SampleLoc,k);
        YQ_avg(j,i)=sum(SampleYQ(index))/k;%get the average channel gain
        %Distributed Kriging 
        [KdB,n_PL,EpsQ]=ChPathLossEsti(SampleYQ(index),distBS); %Estimate the path loss locally
        %Estimate the shadowing correlation locally
        
    end
end