%This file is created by Xu Xiaoli on 14/06/2023
%It generate the channel gain map according to the given channel model. 

clc;
clear;
close all;

%=========Some important parameters=====================
n_PL=2.2;
KdB=-80;
sigma_sq=2;
sigma=sqrt(sigma_sq);
alpha=8;
beta=30;
theta=[KdB; n_PL];

%======MeasureLocations=====
D=500;
d=6;
x_loc=-D/2:d:D/2;
y_loc=-D/2:d:D/2;

[X,Y]=meshgrid(x_loc,y_loc);
[row,col]=size(X);
MeasureLoc=[reshape(X,1,row*col);reshape(Y,1,row*col)];

%Remove the sample at the BS
flag=(MeasureLoc(1,:)==0)+(MeasureLoc(2,:)==0);
%indx=find(flag==2);
MeasureLoc(:,flag==2)=[];

numSamples=length(MeasureLoc(1,:));%the total number of sample points
distance=sqrt(MeasureLoc(1,:).^2+MeasureLoc(2,:).^2); %BS at origin
Hq_all=[ones(numSamples,1),-10*log10(distance')];
distMatrix=sqrt((repmat(MeasureLoc(1,:),numSamples,1)-repmat(MeasureLoc(1,:)',1,numSamples)).^2....
    +(repmat(MeasureLoc(2,:),numSamples,1)-repmat(MeasureLoc(2,:)',1,numSamples)).^2);
RQ_all=alpha*exp(-distMatrix/beta);

YQ=Hq_all*theta+(mvnrnd(zeros(1,numSamples),RQ_all,1))'...
    +(mvnrnd(zeros(1,numSamples),sigma^2*eye(numSamples),1))';
figure;
tri = delaunay(MeasureLoc(1,:)',MeasureLoc(2,:)');
trimesh(tri,MeasureLoc(1,:)',MeasureLoc(2,:)',YQ);


ChannelGainMatrix=reshape(YQ,length(x_loc),length(y_loc));
%==========View the LoS
s=surf(X,Y,ChannelGainMatrix');
s.EdgeColor = 'none';
colorbar;
hold on;
axis equal;

view(2)