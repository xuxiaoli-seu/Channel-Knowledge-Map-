%This file is created by Xu Xiaoli on 02/11/2022
%It simulate the parameter estimation error for PPP and grid distribution,
%and compare it with the analytical results. 
% 
% %========Channel Parameter==================
% alpha=5; %shadowing power
% beta=1; %shadowing correlation distance
% n_PL=3; %path loss exponent
% sigma_sq=4; %multipath variance
% KdB=1; %reference received power
% RegCenter=[200,200]';
% density=25;
% isGrid=1;
% areaSideLength_vec=250;

function [nPLEstiError_simu,KdBEstiError_simu,nPLEstiError_ana,KdBEstiError_ana]=getParaEstiError2(alpha,beta,sigma_sq,RegCenter,isGrid,density,areaSideLength_vec)

%============Region Center with respect to BS=======
%RegCenter=rand(2,1)*D;
dc=sqrt(RegCenter(1)^2+RegCenter(2)^2);
%N_vec=(20:10:300); %number of points

nPLEstiError_simu=zeros(1,length(areaSideLength_vec));
KdBEstiError_simu=zeros(1,length(areaSideLength_vec));
nPLEstiError_ana=zeros(1,length(areaSideLength_vec));
KdBEstiError_ana=zeros(1,length(areaSideLength_vec));
for i=1:length(areaSideLength_vec)
    %N=N_vec(i);
    areaSideLength=areaSideLength_vec(i);
    if isGrid
        d=density;
        coef=max(1,pi*beta^2/d^2);  %expected number of points within shadowing range
        D=areaSideLength;
        MeasureLoc=RegCenter+getMeasureLoc_Grid(D,d);
        N=length(MeasureLoc);
    else
        lambda=density; %get the same density for the Grid distribution
        coef=max(1,pi*beta^2*density); %expected number of points within shadowing range
        D=areaSideLength;
        N=poissrnd( lambda*D^2 );
        MeasureLoc=RegCenter+(rand(2,N)-0.5)*D;       
    end
    %============Simulation Results==========
    distBS=sqrt(MeasureLoc(1,:).^2+MeasureLoc(2,:).^2);
    k=length(distBS);
    distMatrix=sqrt((repmat(MeasureLoc(1,:),k,1)-repmat(MeasureLoc(1,:)',1,k)).^2+(repmat(MeasureLoc(2,:),k,1)-repmat(MeasureLoc(2,:)',1,k)).^2);
    Hq=[ones(k,1),-10*log10(distBS')];
    RQ=alpha*exp(-distMatrix/beta)+sigma_sq*eye(k);
    Tmp=inv(Hq'*Hq);
    CovMat=Tmp*Hq'*RQ*Hq*Tmp;
    nPLEstiError_simu(i)=CovMat(2,2);
    KdBEstiError_simu(i)=CovMat(1,1);
    %=============Analytical Results========
    dmin=max(dc-D/2,0);
    dmax=dc+D/2;
    part1=(dmax*(10*log10(dmax))^2-dmin*(10*log10(dmin))^2)/(dmax-dmin);
    part2=(200*dmax*log10(dmax)-200*dmin*log10(dmin))/(dmax-dmin)/log(10);
    Ki=part1-part2+200/(log(10))^2;
    distVar=100/(log(10))^2-100*dmax*dmin*(log10(dmax/dmin))^2/(dmax-dmin)^2;
    %======
    
    %=======
    nPLEstiError_ana(i)=(alpha+sigma_sq/coef)/(N/coef*distVar);
    KdBEstiError_ana(i)=(alpha+sigma_sq/coef)*Ki/(N/coef*distVar);    
end
% 
% figure;
% plot(N_vec,nPLEstiError_simu,'bo','MarkerFaceColor','b');
% hold on;
% plot(N_vec,nPLEstiError,'b-');
% grid on;
% ylim([0,0.8]);
% xlabel('number of measurements');
% ylabel('Estimation error');
% legend('Simulation','Analysis');
% %title('Estimation Error of the path loss exponent');
% figure;
% plot(N_vec,KdBEstiError_simu,'rs','MarkerFaceColor','r');
% hold on;
% plot(N_vec,KdBEstiError,'r-');
% xlabel('number of measurements');
% ylabel('Estimation error');
% ylim([0,100]);
% legend('Simulation','Analysis');
% title('Estimation Error of K_{dB}');



