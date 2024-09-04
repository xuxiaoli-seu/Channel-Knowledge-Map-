%%This file is created by Xu Xiaoli on 21/02/2023
%IT returns the LS estimation of the Path loss exponent from the input data
%using Least Square Estimation

%It returns the path loss exponent and the residual channel gain (caused by
%shadowing and multipath fading)

function [KdB,n_PL,EpsQ]=ChPathLossEsti(YQ,distBS)

%distBS=sqrt(MeasureLoc(1,:).^2+MeasureLoc(2,:).^2);
k=length(distBS);
Hq=[ones(k,1),-10*log10(distBS')];
tmp=Hq/(Hq'*Hq);
theta_est=YQ*tmp;
KdB=theta_est(1);
n_PL=theta_est(2);
EpsQ=YQ-(Hq*theta_est')';

% if n_PL<0
%     n_PL=0;
%     KdB=sum(YQ)/k; %just take the average
%     EpsQ=YQ-KdB;
% else
%     EpsQ=YQ-(Hq*theta_est')';
% end




