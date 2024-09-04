%This file is created by Xu Xiaoli on 19/10/2023
%It estimates all the channel parameters based on a set of measurements

function [theta,alpha,beta]=estiAllPara(MeasureLoc,YQ)

numSamples=length(YQ); %the number of samples
distance=sqrt(MeasureLoc(1,:).^2+MeasureLoc(2,:).^2);
Hq_sample=[ones(numSamples,1),-10*log10(distance')];
tmp=Hq_sample/(Hq_sample'*Hq_sample);
theta=YQ*tmp;

%Correct the model
if theta(2)<0 % negative path loss exponent
    theta(2)=0;
    theta(1)=sum(YQ)/numSamples;
end

