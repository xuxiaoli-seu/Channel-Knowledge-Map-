%This file is created by Xu Xiaoli on 18/10/2023
%It view the AMSE with estimated paramters when the correlation distance is
%very small

%No need to use the neighboring measurements
clc;
clear;
close all;

d=20;
numSamples_vec=20:50:500;
iter=100;
AMSE_knownPara=zeros(iter,length(numSamples_vec));
AMSE_simu=zeros(iter,length(numSamples_vec));

%===============================================
%Generate the data
%===============================================
%Specify the channel parameters 
alpha=8; %shadowing power
beta=1; %shadowing correlation distance
n_PL=2.2; %path loss exponent
sigma_sq=2; %multipath variance
KdB=-80; %reference received power
theta=[KdB; n_PL];
AMSE_lowBeta_ana=(alpha+sigma_sq)*(1+2./numSamples_vec);
%RegionCenter=[200;200];

for j=1:iter
    j
    RegionCenter=500*rand(2,1);
    for i=1:length(numSamples_vec)
        numSamples_lb=numSamples_vec(i);
        D=sqrt(numSamples_lb)*d;
        %MeasureLoc=getMeasureLoc_Grid_GivenCenter(D,d,RegionCenter);
        density=numSamples_lb/D^2;
        MeasureLoc=getMeasureLoc_PPP_GivenCenter(D,density,RegionCenter);
        
        numTarget=2000; 
        TargetLocation=RegionCenter+[(rand(1,numTarget)-1/2)*D;(rand(1,numTarget)-1/2)*D]; %get the random locations within the region
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

        tmp=Hq_sample/(Hq_sample'*Hq_sample);
        theta_est=YQ_sample'*tmp;
    %     KdB_esti=theta_est(1);
    %     n_PL_esti=theta_est(2);
        %YQ_shaFad=YQ_sample-(Hq_sample*theta_est');
        %[alpha_est1, beta_est1,sigma_est1]=EstShadowPara_LS(YQ_shaFad,distanceMatrixSample);
        %Estimate the paramters
        Hq_target=Hq_all(numSamples+1:end,:);
        YQ_target_estiPara=Hq_target*theta_est';
        YQ_target_knownPara=Hq_target*theta;
        AMSE_simu(j,i)=sum((YQ_target_estiPara-YQ_target).^2)/numTarget;
        AMSE_knownPara(j,i)=sum((YQ_target_knownPara-YQ_target).^2)/numTarget;
    end
end

figure;
plot(numSamples_vec,AMSE_lowBeta_ana,'k--');
hold on;
plot(numSamples_vec,sum(AMSE_simu,1)/iter,'ko-');
plot(numSamples_vec,sum(AMSE_knownPara,1)/iter,'rs-');


% plot(k_vec,sum(AMSE_knownPara3)/iter,'b--');
% plot(k_vec,sum(AMSE_simu3)/iter,'bo-');
% plot(k_vec,sum(AMSE_lowBeta_ana3)/iter*ones(1,length(k_vec)),'b*-');


%plot(k_vec,AMSE_simu2,'rs-');
hold off;
xlabel('$k$','interpreter','latex');
ylabel('AMSE');
legend('Analysis','Estimated Para','Known Para');
save LowBetaAMSE_PPP.mat;