%This file is created by Xu Xiaoli on 08/09/2023
%It estimates the shadowing power and multipath fading variance of the
%channel--based on the simulated data

%Two methods are considered: 
%1. LS estimation 
%2. Semivariogram Fitting

clc;
clear;
close all;

%===============================================
%Generate the data
%===============================================
%Specify the channel parameters 
alpha=8; %shadowing power
beta=30; %shadowing correlation distance
n_PL=2.2; %path loss exponent
sigma_sq=2; %multipath variance
KdB=-80; %reference received power
theta=[KdB; n_PL];
% The Measurements and Test Data set
numSamples_lb=30; % the minimum number of samples in the CGM
RegionCenter=[100;100];
d=25; %sample separation
D=sqrt(numSamples_lb)*d; % determine the region size to ensure there are at least 10 times the samples 
% The locations of the Measurement points and Test points
MeasureLoc=getMeasureLoc_Grid_GivenCenter(D,d,RegionCenter);
numTarget=5000; %the number of test locations
TargetLocation=RegionCenter+[(rand(1,numTarget)-1/2)*D;(rand(1,numTarget)-1/2)*D]; %get the random locations within the region

%View the sample locations
% figure;
% scatter(MeasureLoc(1,:),MeasureLoc(2,:),'o','MarkerFaceColor','b');
% hold on;
% scatter(TargetLocation(1,:),TargetLocation(2,:),'o','MarkerFaceColor','r');

AllPointLoc=[MeasureLoc,TargetLocation];
numSamples=length(MeasureLoc(1,:));%the total number of sample points
numTarget=length(TargetLocation);
Totalnum=numSamples+numTarget;
distance=sqrt(AllPointLoc(1,:).^2+AllPointLoc(2,:).^2); %BS at origin
Hq_all=[ones(Totalnum,1),-10*log10(distance')];
Hq_sample=[ones(numSamples,1),-10*log10(distance(1:numSamples)')];

distMatrixAll=sqrt((repmat(AllPointLoc(1,:),Totalnum,1)-repmat(AllPointLoc(1,:)',1,Totalnum)).^2....
    +(repmat(AllPointLoc(2,:),Totalnum,1)-repmat(AllPointLoc(2,:)',1,Totalnum)).^2); %the distance matrix for both samples and targets
distanceMatrixSample=sqrt((repmat(MeasureLoc(1,:),numSamples,1)-repmat(MeasureLoc(1,:)',1,numSamples)).^2....
    +(repmat(MeasureLoc(2,:),numSamples,1)-repmat(MeasureLoc(2,:)',1,numSamples)).^2); %the distance matrix only for the samples

RQ_all=alpha*exp(-distMatrixAll/beta);

YQ_all=Hq_all*theta+(mvnrnd(zeros(1,Totalnum),RQ_all,1))'...
    +(mvnrnd(zeros(1,Totalnum),sigma_sq*eye(Totalnum),1))';
YQ_sample=YQ_all(1:numSamples);
YQ_target=YQ_all(numSamples+1:end);

% View the obtained channel gain
% figure;
% tri = delaunay(MeasureLoc(1,:)',MeasureLoc(2,:)');
% trimesh(tri,MeasureLoc(1,:)',MeasureLoc(2,:)',YQ);

%==================================================
%Estimate the path loss parameters based on the measurements
%==================================================
tmp=Hq_sample/(Hq_sample'*Hq_sample);
theta_est=YQ_sample'*tmp;
KdB_esti=theta_est(1);
n_PL_esti=theta_est(2);

%===================================================
%Estimate the shadowing power and mutlipath variance
% Method 1: Least Square
% Method 2: Variogram Fitting. 
%===================================================
YQ_shaFad=YQ_sample-(Hq_sample*theta_est'); %get the shadowing and fading part

[alpha_est1, beta_est1,sigma_est1]=EstShadowPara_LS(YQ_shaFad,distanceMatrixSample);
%[alpha_est2, beta_est2,sigma_est2]=EstShadowPara_VariFit(YQ_shaFad,MeasureLoc);

%======================================================
% Simulate the AMSE based on the estimated parameters
%======================================================

%Test the estimation when k=0, i.e., only the path loss is counted
% YQ_esti=-10*log10(distance(numSamples+1:end))*n_PL+KdB; %distance to the BS
% MSE=sum((YQ_esti'-YQ_all(numSamples+1:end)).^2)/numTarget
% 
% 
% return;


k_vec=0:10;
AMSE_knownPara=zeros(1,length(k_vec));
AMSE_simu1=zeros(1,length(k_vec));
AMSE_simu2=zeros(1,length(k_vec));
for i=1:length(k_vec)
    k=k_vec(i);
    AMSE_knownPara(i)=AMSEsimu_estiPara(MeasureLoc,TargetLocation,YQ_all,k,n_PL,KdB,alpha,beta);
    AMSE_simu1(i)=AMSEsimu_estiPara(MeasureLoc,TargetLocation,YQ_all,k,n_PL_esti,KdB_esti,alpha_est1,beta_est1);
    %AMSE_simu2(i)=AMSEsimu_estiPara(MeasureLoc,TargetLocation,YQ_all,k,n_PL_esti,KdB_esti,alpha_est2,beta_est2);
end

figure;
plot(k_vec,AMSE_knownPara,'ko--');
hold on;
plot(k_vec,AMSE_simu1,'bo-');
%plot(k_vec,AMSE_simu2,'rs-');
hold off;
xlabel('$k$','interpreter','latex');
ylabel('AMSE');
legend('Known Para','Estimated Para');

%legend('Known Para','Estimated Para (LS)', 'Estimated Para (Vario Fitting)');



