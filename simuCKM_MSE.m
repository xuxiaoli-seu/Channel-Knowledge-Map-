%This file is created by Xu Xiaoli on 13/06/2022
%It simulate the average normalized MSE by assuming Grid/PPP distribution
%of the measurement positions. 
%The esimtaiton of randomly selected location is assumed to be the MMSE
%estimator 
%Input:
%D: the sidelength of the area
%density: the number of meansurements in unit area (1m^2)


%k: k-nearest neighbore used for estimation
%flag: if flag==0, grid sample is used; if flag==1, PPP sample is used

function AMSE=simuCKM_MSE(D,k,MeasureLoc)
%===========Envrionment Loc==========================
%BS=[0,0]; %BS locate at the center of the region
%D=1000; %the sidelength of the area;

%========Channel Parameter==================
global alpha; %shadowing power
global beta; %shadowing correlation distance
global n_PL; %path loss exponent
global sigma; %multipath variance
global KdB; %reference received power

theta=[KdB; n_PL];

%=============Simulate the measurement results===
numSamples=length(MeasureLoc(1,:));%the total number of sample points
distance=sqrt(MeasureLoc(1,:).^2+MeasureLoc(2,:).^2); %BS at origin
Hq_all=[ones(numSamples,1),-10*log10(distance')];
distMatrix=sqrt((repmat(MeasureLoc(1,:),numSamples,1)-repmat(MeasureLoc(1,:)',1,numSamples)).^2....
    +(repmat(MeasureLoc(2,:),numSamples,1)-repmat(MeasureLoc(2,:)',1,numSamples)).^2);
RQ_all=alpha*exp(-distMatrix/beta);

YQ=Hq_all*theta+(mvnrnd(zeros(1,numSamples),RQ_all,1))'...
    +(mvnrnd(zeros(1,numSamples),sigma^2*eye(numSamples),1))';
% figure;
% tri = delaunay(MeasureLoc(1,:)',MeasureLoc(2,:)');
% trimesh(tri,MeasureLoc(1,:)',MeasureLoc(2,:)',YQ);

%===========================================
%k_vec=1:10; %number of measurements used for estimate the channel at q
SampleSize=1000;
EstiPower=zeros(1,SampleSize);
MSE=zeros(1,SampleSize);

for i=1:SampleSize
    q_Loc=[(rand-1/2)*D;(rand-1/2)*D];
    [distk,distMatrix,distBS,idx]=findkNearest(q_Loc,MeasureLoc,k);
    Hq=[ones(k,1),-10*log10(distBS')];%2*k,distance from the selected measurements to BS
    hq=[1,-10*log10(norm(q_Loc))];

    PhiQ=alpha*exp(-distk'/beta);
    RQ=alpha*exp(-distMatrix/beta)+sigma^2*eye(k);
    tmp=PhiQ'/RQ;
    EstiPower(i)=hq*theta+tmp*(YQ(idx)-Hq*theta);
    MSE(i)=alpha+sigma^2-tmp*PhiQ;
end


AMSE=sum(MSE,2)/length(MSE); %normoalized mean sequre error
% figure;
% plot(k_vec,ANMSE,'bo-','MarkerFaceColor','b');
% xlabel('Number of Neighbors, k');
% ylabel('average MSE');

