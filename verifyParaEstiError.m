%This file is created by Xu Xiaoli on 02/11/2022
%It simulate the parameter estimation error for PPP and grid distribution,
%and compare it with the analytical results. 
clc;
clear;
close all;

%========Channel Parameter==================
alpha=8; %shadowing power
beta=30; %shadowing correlation distance
n_PL=2.2; %path loss exponent
sigma_sq=4; %multipath variance
KdB=-80; %reference received power
theta=[KdB; n_PL];

%=============CKM density and size============
BS=[0,0]; %BS locate at the center of the region
D=500; %the sidelength of the area;
lambda=1e-3;
d=sqrt(1/lambda);%the average separateion in grid distribution
d=10;
isGrid=1; %if 1, it is grid, otherwise it is PPP

%============Region Center with respect to BS=======
%RegCenter=rand(2,1)*D;
RegCenter=[100,100]';
dc=sqrt(RegCenter(1)^2+RegCenter(2)^2);
N_vec=(20:10:300); %number of points
nPLEstiError_simu=zeros(1,length(N_vec));
KdBEstiError_simu=zeros(1,length(N_vec));

nPLEstiError=zeros(1,length(N_vec));
KdBEstiError=zeros(1,length(N_vec));
for i=1:length(N_vec)
    N=N_vec(i);
    
    if isGrid
        D=sqrt(N)*(d-1);
        MeasureLoc=RegCenter+getMeasureLoc_Grid(D,d);
        
    else
        D=sqrt(N/lambda);
        MeasureLoc=RegCenter+(rand(2,N)-0.5)*D;
        
    end
    %============Simulation Results==========
    distBS=sqrt(MeasureLoc(1,:).^2+MeasureLoc(2,:).^2);
    dmin=max(dc-D/2,0);
    dmax=dc+D/2;
    k=length(distBS);
    distMatrix=sqrt((repmat(MeasureLoc(1,:),k,1)-repmat(MeasureLoc(1,:)',1,k)).^2+(repmat(MeasureLoc(2,:),k,1)-repmat(MeasureLoc(2,:)',1,k)).^2);
    Hq=[ones(k,1),-10*log10(distBS')];
    RQ=alpha*exp(-distMatrix/beta)+sigma_sq*eye(k);
    Tmp=inv(Hq'*Hq);
    CovMat=Tmp*Hq'*RQ*Hq*Tmp;
    nPLEstiError_simu(i)=CovMat(2,2)
    KdBEstiError_simu(i)=CovMat(1,1);
    %===========Analytical Results============
    distMean=(10*dmax*log10(dmax)-10*dmin*log10(dmin))/(dmax-dmin)-10/log(10);
    part1=(dmax*(10*log10(dmax))^2-dmin*(10*log10(dmin))^2)/(dmax-dmin);
    part2=(200*dmax*log10(dmax)-200*dmin*log10(dmin))/(dmax-dmin)/log(10);
    Ki=part1-part2+200/(log(10))^2;
    distVar=100/(log(10))^2-100*dmax*dmin*(log10(dmax/dmin))^2/(dmax-dmin)^2;

    nPLEstiError(i)=(alpha+sigma_sq)/(N*distVar);
    KdBEstiError(i)=(alpha+sigma_sq)*Ki/(N*distVar);
end

figure;
plot(N_vec,nPLEstiError_simu,'bo','MarkerFaceColor','b');
hold on;
plot(N_vec,nPLEstiError,'b-');
grid on;
ylim([0,0.8]);
xlabel('number of measurements');
ylabel('Estimation error');
legend('Simulation','Analysis');
%title('Estimation Error of the path loss exponent');
figure;
plot(N_vec,KdBEstiError_simu,'rs','MarkerFaceColor','r');
hold on;
plot(N_vec,KdBEstiError,'r-');
xlabel('number of measurements');
ylabel('Estimation error');
ylim([0,100]);
legend('Simulation','Analysis');
title('Estimation Error of K_{dB}');



