%This file is created by Xu Xiaoli on 13/06/2022
%It returns the measurement positions for Grid distribution
% 
% D=1000; %the sidelength of the area;
% d=10; %the average separateion in grid distribution

function MeasureLoc=getMeasureLoc_Grid_GivenCenter(D,d,RegionCenter)

x_loc=-D/2:d:D/2;
y_loc=-D/2:d:D/2;

[X,Y]=meshgrid(x_loc,y_loc);
[row,col]=size(X);
MeasureLoc=[reshape(X,1,row*col);reshape(Y,1,row*col)]+RegionCenter;

%Remove the sample at the BS
flag=(MeasureLoc(1,:)==0)+(MeasureLoc(2,:)==0);
%indx=find(flag==2);
MeasureLoc(:,flag==2)=[];
    
%view the distribution of measurement loc;
% scatter(MeasureLoc(1,:),MeasureLoc(2,:),'o','MarkerFaceColor','b');
% hold on;
% plot(0,0,'rp','MarkerFaceColor','r');
% hold off;