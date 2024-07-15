%This file is created by Xu Xiaoli on 22/06/2022
%It return the simulated MSE.

clc;
clear;
close all;
global D;
global alpha
global beta
global n_PL
global sigma
global KdB

D=1000;

%MeasureLoc=getMeasureLoc_PPP(D,d);

d=25;
MeasureLoc=getMeasureLoc_Grid(D,d);

alpha=8; %shadowing power
beta=30; %shadowing correlation distance
n_PL=2.2; %path loss exponent
sigma=sqrt(2); %multipath variance
KdB=-80; %reference received power

iter=200;
k_vec=0:10;
AvgMSE_simu=zeros(length(k_vec),iter);
AvgMSE_ana=zeros(length(k_vec),iter);

for i=1:iter
    for j=1:length(k_vec)
        k=k_vec(j)
        [AvgMSE_simu(j,i), AvgMSE_ana(j,i)]=getMSEsimu(MeasureLoc,100,k);
    end
end

save Grid_d25.mat;

plot(k_vec,sum(AvgMSE_simu,2)/iter,'ro');
hold on;
plot(k_vec,sum(AvgMSE_ana,2)/iter,'b--');
plot(k_vec,(alpha+sigma^2-alpha^2/(alpha+sigma^2))*ones(1,length(k_vec)),'k--');
hold off;
xlabel('k');
ylabel('Average MSE');
legend('simu','ana');
