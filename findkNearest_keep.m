%This file is created by Xu Xiaoli on 04/03/2024
%It returns the distance vector between q_loc and k nearest neighbors
% Only keep the sample distance no greater than the correlation distance to
% ensure that the matrix inverse is feasible.

%distk: is the vector store the distance between the target and the nearest
%k measurements

%distMatrix: k*k matrix store the distance between any two measurement
%matrix
function [distk,distMatrix,distBS,idx]=findkNearest_keep(q_Loc,MeasureLoc,k, beta)

dist_vec=sqrt((MeasureLoc(1,:)-q_Loc(1)).^2+(MeasureLoc(2,:)-q_Loc(2)).^2);

[distk,idx0]=mink(dist_vec,k);
keepIdx=find(distk<beta,1,'last');
idx=idx0(1:keepIdx);
selectedLoc=MeasureLoc(:,idx)
k=length(idx)
%distance between the measurements
distMatrix=sqrt((repmat(selectedLoc(1,:),k,1)-repmat(selectedLoc(1,:)',1,k)).^2+(repmat(selectedLoc(2,:),k,1)-repmat(selectedLoc(2,:)',1,k)).^2);

%distance to the BS
distBS=sqrt(selectedLoc(1,:).^2+selectedLoc(2,:).^2);

