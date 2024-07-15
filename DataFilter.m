%This file is created by Xu Xiaoli on 22/02/2023
%It select the reasonable value with  Nh>30 and k<H/2
% It applies the practical rule for the weighted LS esimation proposed in 
% [1]N. Cressie, "Fitting variogram models by weighted least squares", J.
% Int. Assoc. Math. Geol. vol. 17, pp. 563-586, 1985. 

function [h,Nh,Vario]=DataFilter(h,Nh,Vario) 

%H=min(max(h),100);
indx=find(Nh>=30);
indx2=find(h>200,1,'first');
%k=min(floor(H/2),length(indx));
k=min(indx2,length(indx));
h=h(indx(1:k));
Nh=Nh(indx(1:k));
Vario=Vario(indx(1:k));

