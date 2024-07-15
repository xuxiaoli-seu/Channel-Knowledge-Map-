%This file is created by Xu Xiaoli on 04/03/2024
%It returns the measurement positions for PPP distribution
% 


function MeasureLoc=getMeasureLoc_PPP_GivenCenter(D,numSamples,RegionCenter)

N = numSamples;
x_loc=(rand(1,N)-1/2)*D;
y_loc=(rand(1,N)-1/2)*D;
MeasureLoc=[x_loc;y_loc]+RegionCenter;




%Remove the sample at the BS
% flag=(MeasureLoc(1,:)==0)+(MeasureLoc(2,:)==0);
% %indx=find(flag==2);
% MeasureLoc(:,flag==2)=[];
    
%view the distribution of measurement loc;
% scatter(MeasureLoc(1,:),MeasureLoc(2,:),'o','MarkerFaceColor','b');
% hold on;
% plot(0,0,'rp','MarkerFaceColor','r');
% hold off;