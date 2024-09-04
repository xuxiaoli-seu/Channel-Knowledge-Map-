%This file is created by Xu Xiaoli on 27/09/2022
%It read the map from txt file

clc;
clear;
close all;

TxLoc=[-482.27,-642.89,100];
fileID = fopen('CKM.txt','r');
CKMdata=fscanf(fileID,['%f','%f','%f','%f','%f']);
%[Idx,RxLoc1,RxLoc2,RxLoc3,Gain,LoS]=fscanf(fileID,['%d','%f','%f','%f','%d']);
numPoints=length(CKMdata)/6;
CKMReshape=reshape(CKMdata,[6,numPoints]);
%MeasureLoc=CKMReshape(2:4,:);
ChannelGain=CKMReshape(5,:);
LoSFlag=CKMReshape(6,:);

ChannelGainMatrix=reshape(ChannelGain,401,401);
XRange=-610:1:-210;
YRange=-810:1:-410;
%=====Shift the Data to put the Tx at the origin=========
XLoc=XRange-TxLoc(1);
YLoc=YRange-TxLoc(2);
LoSMatrix=reshape(LoSFlag,401,401);


[X,Y]=meshgrid(XLoc,YLoc);
%==========View the LoS
s=surf(X,Y,ChannelGainMatrix');
s.EdgeColor = 'none';
colorbar;
hold on;

view(2)
xlim([-150,280]);
ylim([-180,250]);
plot(0,0,'rp','MarkerFaceColor','r');

%=========Extract the data for parameter estimation=======================
d=20;
xlength=length(XLoc);
selectedPointIdx=1:d:xlength;
SelectedXLoc=XLoc(selectedPointIdx);
SelectedYLoc=YLoc(selectedPointIdx);
SelectedMeasurement=ChannelGainMatrix(selectedPointIdx,selectedPointIdx)';
SelectedLoS=LoSMatrix(selectedPointIdx,selectedPointIdx)';
[Xs,Ys]=meshgrid(SelectedXLoc,SelectedYLoc);

remainPointIdx=setdiff(1:xlength,selectedPointIdx);
remainXLoc=XLoc(remainPointIdx);
remainYLoc=YLoc(remainPointIdx);
remainMeasurement=ChannelGainMatrix(remainPointIdx,remainPointIdx)';
remainLoS=LoSMatrix(remainPointIdx,remainPointIdx)';
[Xr,Yr]=meshgrid(remainXLoc,remainYLoc);

%Plot the location of the selected measurements
% figure;
% plot(X,Y,'bo','MarkerFaceColor','b','MarkerSize',1);
% xlim([-150,280]);
% ylim([-180,250]);
% figure;
% s=surf(X,Y,remainMeasurement);
% s.EdgeColor = 'none';
% view(2)
% xlim([-150,280]);
% ylim([-180,250]);
% hold on;
% plot(X,Y,'bo','MarkerFaceColor','b','MarkerSize',1);

%============================

[row,col]=size(Xs);
MeasureLoc_all=[reshape(Xs,1,row*col);reshape(Ys,1,row*col)];
YQ_all=reshape(SelectedMeasurement,1,row*col);
LoS_all=reshape(SelectedLoS,1,row*col);

Indx_keep=(YQ_all>-200); %drop the points with recieved power less than -200dB
MeasureLoc=MeasureLoc_all(:,Indx_keep);
YQ=YQ_all(Indx_keep);
LoS=LoS_all(Indx_keep);



[row,col]=size(Xr);
MeasureLoc_remain_all=[reshape(Xr,1,row*col);reshape(Yr,1,row*col)];
YQ_remain_all=reshape(remainMeasurement,1,row*col);
LoS_remain_all=reshape(remainLoS,1,row*col);

Indx_keep=(YQ_remain_all>-200); %drop the points with recieved power less than -200dB
MeasureLoc_remain=MeasureLoc_remain_all(:,Indx_keep);
YQ_remain=YQ_remain_all(Indx_keep);
LoS_remain=LoS_remain_all(Indx_keep);

figure;
plot(MeasureLoc(1,:),MeasureLoc(2,:),'bo','MarkerFaceColor','b','MarkerSize',2);
hold on;
plot(0,0,'rp','MarkerFaceColor','r');
xlim([-140,250]);
ylim([-165, 250]);
%===========Remove the Useless points (with Rx<-200dB)======


save MeasureData2_4GHz_d20.mat;


% figure;
% tri = delaunay(MeasureLoc(1,:)',MeasureLoc(2,:)');
% trimesh(tri,MeasureLoc(1,:)',MeasureLoc(2,:)',YQ);
