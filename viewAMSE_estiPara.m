clc;
clear;
close all;
iter=100;
k_vec=0:1:8;
d=25;
density=1/d^2;
numSamples_vec=50;
AMSE_knownPara1=zeros(iter,length(k_vec));
AMSE_simu1=zeros(iter,length(k_vec));
AMSE_lowBeta_ana1=zeros(iter,1);

AMSE_knownPara2=zeros(iter,length(k_vec));
AMSE_simu2=zeros(iter,length(k_vec));
AMSE_lowBeta_ana2=zeros(iter,1);

AMSE_knownPara3=zeros(iter,length(k_vec));
AMSE_simu3=zeros(iter,length(k_vec));
AMSE_lowBeta_ana3=zeros(iter,1);

for i=1:iter
    i
    RegionCenter=500*rand(2,1);
    [AMSE_knownPara1(i,:),AMSE_simu1(i,:),AMSE_lowBeta_ana1(i)]=getAMSEwithEstiPara(RegionCenter,d,k_vec,numSamples_vec(1));
end
MSEana=MSEPPPana(density,k_vec,8,30,sqrt(2));
%MSEana=MSEGridana(d,k_vec,8,30,sqrt(2));

figure;
plot(k_vec,sum(AMSE_knownPara1,1)/iter,'rs-');
hold on;
plot(k_vec,sum(AMSE_simu1,1)/iter,'ko--');
plot(k_vec,MSEana,'r--');
%plot(k_vec,sum(AMSE_lowBeta_ana1)/iter*ones(1,length(k_vec)),'k*-');
% 
% plot(k_vec,sum(AMSE_knownPara2,1)/iter,'rs-');
% plot(k_vec,sum(AMSE_simu2,1)/iter,'ko-');
%plot(k_vec,sum(AMSE_lowBeta_ana2)/iter*ones(1,length(k_vec)),'r*-');

% 
% plot(k_vec,sum(AMSE_knownPara3,1)/iter,'b--');
% plot(k_vec,sum(AMSE_simu3,1)/iter,'bo-');
% %plot(k_vec,sum(AMSE_lowBeta_ana3)/iter*ones(1,length(k_vec)),'b*-');


%plot(k_vec,AMSE_simu2,'rs-');

xlabel('$k$','interpreter','latex');
ylabel('AMSE');
legend('Known Para.','Estimated Para.','Analysis');

% clear;
% load AMSEestiParaLargeBeta_d25.mat

tmp1=sum(AMSE_knownPara1,1)/iter;
tmp2=sum(AMSE_simu1,1)/iter;
MSE=MSEGridana(d,k_vec,8,30,sqrt(2));
% create a new pair of axes inside current figure
axes('position',[.65 .175 .25 .25])
box on % put box around new pair of axes
indexOfInterest = 1:4; % range of t near perturbation
plot(k_vec(indexOfInterest),tmp1(indexOfInterest),'rs-','LineWidth',1) % plot on new axes
hold on;
plot(k_vec(indexOfInterest),tmp2(indexOfInterest),'ko-','LineWidth',1) % plot on new axes
plot(k_vec(indexOfInterest),MSE(indexOfInterest),'r--','LineWidth',1) % plot on new axes
hold off;

axis tight

save AMSEestiParaBeta30_PPP_d25_N100.mat
