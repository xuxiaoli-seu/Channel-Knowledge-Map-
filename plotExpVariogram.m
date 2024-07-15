%This file is created by Xu Xiaoli on 21/02/2023
%It plots the experimental variogram of the data input 

%The input data is assumed to be intrinsic stationary

%MeasureLoc: the set of locations for the data
%Data: the measured data at the locations




function [h,Nh,Varioh,RobustVario]=plotExpVariogram(MeasureLoc,Data)

numSamples=length(MeasureLoc(1,:));%the total number of sample points
distMatrix=sqrt((repmat(MeasureLoc(1,:),numSamples,1)-repmat(MeasureLoc(1,:)',1,numSamples)).^2....
    +(repmat(MeasureLoc(2,:),numSamples,1)-repmat(MeasureLoc(2,:)',1,numSamples)).^2);
distRound=round(distMatrix);
distSample=reshape(distRound,1,numSamples*numSamples);
h=unique(distSample);
h=h(2:end);
Nh=zeros(1,length(h));
TotalVari=zeros(1,length(h));
TotalRobust=zeros(1,length(h));

for i=1:numSamples
    for j=i+1:numSamples
        hIdx=find(h==distRound(i,j));
        Nh(hIdx)=Nh(hIdx)+1;
        TotalVari(hIdx)=TotalVari(hIdx)+(Data(i)-Data(j))^2;
        TotalRobust(hIdx)=TotalRobust(hIdx)+sqrt(abs(Data(i)-Data(j)));
    end
end
Varioh=TotalVari./Nh/2;
RobustVario=((TotalRobust./Nh).^4)./(0.457+0.494./Nh)/2;

% figure;
% plot(h,Varioh,'bo');





