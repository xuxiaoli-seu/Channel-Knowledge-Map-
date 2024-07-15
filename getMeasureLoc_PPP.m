%This file is created by Xu Xiaoli on 13/06/2022
%It returns the measurement positions for Grid distribution
% 
% D=1000; %the sidelength of the area;
% d=10; %the average separateion in grid distribution

function MeasureLoc=getMeasureLoc_PPP(D,density)

N = poissrnd( density*D^2 ) ;
x_loc=(rand(1,N)-1/2)*D;
y_loc=(rand(1,N)-1/2)*D;

MeasureLoc=[x_loc;y_loc];

%Remove the sample at the BS
flag=(MeasureLoc(1,:)==0)+(MeasureLoc(2,:)==0);
%indx=find(flag==2);
MeasureLoc(:,flag==2)=[];
    
%view the distribution of measurement loc;
figure;
scatter(MeasureLoc(1,:),MeasureLoc(2,:),'o','SizeData',5,'MarkerFaceColor','b');
hold on;
plot(0,0,'rp','MarkerFaceColor','r');
hold off;