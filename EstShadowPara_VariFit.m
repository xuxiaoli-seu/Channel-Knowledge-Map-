%This file is created by Xu Xiaoli on 07/10/2023
%It estimates the shadowing parameters for based on the variogram fitting 

%alpha: shadowing variance
%beta: correlation distance
%sigma: multipath fading variance

%YQ_shaFad: the channel gain substract the path loss, equals to
%shadowing+fading

function [alpha, beta, sigma_sq]=EstShadowPara_VariFit(YQ_shaFad,MeasureLoc)

%k=length(YQ_shaFad); % the number of measurements

%Chi_esti=(YQ_shaFad'*YQ_shaFad)/k; %the total variance, equals alpha+sigma^2

[h,Nh,Varioh,RobustVario]=plotExpVariogram(MeasureLoc,YQ_shaFad);
[h,Nh,Vario]=DataFilter(h,Nh,Varioh);
[VarioPara,~]=WeightedLSFit(h,Nh,Vario,3);%use the exponential model
alpha=VarioPara(2);
sigma_sq=VarioPara(1);
beta=VarioPara(3);

