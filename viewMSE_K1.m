%This file is created by Xu Xiaoli on 13/06/2022
%It simulate the average normalized MSE by assuming Grid/PPP distribution
%of the measurement positions. 
%The esimtaiton of randomly selected location is assumed to be the MMSE



clc;
clear;
close all;
global alpha
global beta
global n_PL
global sigma
global KdB


alpha=8; %shadowing power
beta=30; %shadowing correlation distance
n_PL=2.2; %path loss exponent
sigma=sqrt(2); %multipath variance
KdB=-80; %reference received power
density_vec=(1:5:40)*1e-4;
%density_vec=[1e-2 1.6e-3];

iter=10;
MSE_vec_PPP=zeros(iter,length(density_vec));
MSE_vec_grid=zeros(iter,length(density_vec));
k=1;
D=1000;
d_vec=sqrt(1./density_vec);
for j=1:iter
    for i=1:length(density_vec)
        density=density_vec(i);
        d=d_vec(i);

        MeasureLoc=getMeasureLoc_PPP(D,density);
        MSE_vec_PPP(j,i)=simuCKM_MSE(D, k, MeasureLoc);

        MeasureLoc=getMeasureLoc_Grid(D,d);
        MSE_vec_grid(j,i)=simuCKM_MSE(D, k, MeasureLoc);
    end
end
MSE_vec_PPP_avg=sum(MSE_vec_PPP)/iter;
MSE_vec_grid_avg=sum(MSE_vec_grid)/iter;

coef=pi*density_vec*beta^2;
AMSE_ana_PPP=alpha+sigma^2-alpha^2/(alpha+sigma^2)*...
    (1-1./(beta*sqrt(density_vec)).*exp(1./coef).*(1-erf(1./sqrt(coef))));


I1=1-(1+sqrt(2)*d_vec/beta).*exp(-sqrt(2)*d_vec/beta);
%Verify the second integration
I2_simu=zeros(1,length(d_vec));
for i=1:length(d_vec)
    d=d_vec(i);
    x=d/2:0.01:sqrt(2)*d/2;
    y=(8*x/d^2).*exp(-2*x/beta).*acos(0.5*d./x);
    I2_simu(i)=sum(y)/100;
end
I2_ana=0.5708-0.7158*d_vec/beta;
AMSE_ana_Grid=alpha+sigma^2-alpha^2/(alpha+sigma^2)*((pi*beta^2./(2*d_vec.^2)).*I1-I2_simu);

a=d_vec/beta;
I2_approx=2*exp(-a).*(0.2854-0.0725*a+0.0108*a.^2);
AMSE_approx_Grid=alpha+sigma^2-alpha^2/(alpha+sigma^2)*((pi*beta^2./(2*d_vec.^2)).*I1-I2_approx);

figure;
plot(d_vec, AMSE_ana_Grid,'r--');
hold on;
plot(d_vec,MSE_vec_grid_avg,'rs');
plot(d_vec,AMSE_approx_Grid,'k*');
hold off;
xlabel('Separation distance d');
ylabel('Average MSE');
legend('Grid-ana','Grid-simu','Grid-approx (closed-form)');

figure;
plot(d_vec, AMSE_ana_Grid,'r--');
hold on;
plot(d_vec,MSE_vec_grid_avg,'rs','MarkerFaceColor','r');
plot(density_vec,AMSE_ana_PPP,'b--');
plot(density_vec,MSE_vec_PPP_avg,'bo','MarkerFaceColor','b');

hold off;
xlabel('Density (per m^2)');
ylabel('Average MSE');
legend('Grid-ana','Grid-simu','PPP-ana','PPP-simu');
grid on;


