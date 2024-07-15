%This file is created by Xu Xiaoli on 13/06/2022
%It returns the distance vector between q_loc and k nearest neighbors

%distk: is the vector store the distance between the target and the nearest
%k measurements

%distMatrix: k*k matrix store the distance between any two measurement
%matrix
function [distk,distMatrix,distBS,idx]=findkNearest(q_Loc,MeasureLoc,k)

dist_vec=sqrt((MeasureLoc(1,:)-q_Loc(1)).^2+(MeasureLoc(2,:)-q_Loc(2)).^2);

[distk,idx]=mink(dist_vec,k);

selectedLoc=MeasureLoc(:,idx);
%distance between the measurements
distMatrix=sqrt((repmat(selectedLoc(1,:),k,1)-repmat(selectedLoc(1,:)',1,k)).^2+(repmat(selectedLoc(2,:),k,1)-repmat(selectedLoc(2,:)',1,k)).^2);

%distance to the BS
distBS=sqrt(selectedLoc(1,:).^2+selectedLoc(2,:).^2);

