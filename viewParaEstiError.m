%This file is created by Xu Xiaoli on 15/11/2022
%It compares the performance of the parameter estimation with the
%analytical results

clc;
clear;
close all;


%========Channel Parameter==================
alpha=8; %shadowing power
beta=30; %shadowing correlation distance
%n_PL=2.2; %path loss exponent
sigma_sq=2; %multipath variance
%theta=[KdB; n_PL];

%=============CKM density and size============
BS=[0,0]; %BS locate at the center of the region
D=1000; %the sidelength of the area;
%lambda=1/d^2; %density for PPP samples
%isGrid=1; %if 1, it is grid, otherwise it is PPP

%=====================================
RegCenter=[200,200]'; %the center region
areaSideLength_vec=150:50:500;

%N_vec=30:10:500;
d=5;
[nPLEstiError_simu,KdBEstiError_simu,nPLEstiError_ana,KdBEstiError_ana]=getParaEstiError2(alpha,beta,sigma_sq,RegCenter,1,d,areaSideLength_vec);
[nPLEstiError_simu2,KdBEstiError_simu2,nPLEstiError_ana2,KdBEstiError_ana2]=getParaEstiError2(alpha,beta,sigma_sq,RegCenter,0,1/d^2,areaSideLength_vec);
d=10;
[nPLEstiError_simu3,KdBEstiError_simu3,nPLEstiError_ana3,KdBEstiError_ana3]=getParaEstiError2(alpha,beta,sigma_sq,RegCenter,1,d,areaSideLength_vec);
[nPLEstiError_simu4,KdBEstiError_simu4,nPLEstiError_ana4,KdBEstiError_ana4]=getParaEstiError2(alpha,beta,sigma_sq,RegCenter,0,1/d^2,areaSideLength_vec);

figure;
semilogy(areaSideLength_vec,nPLEstiError_simu,'bo');
hold on;
semilogy(areaSideLength_vec,nPLEstiError_ana,'b--','LineWidth',1);
semilogy(areaSideLength_vec,nPLEstiError_simu2,'rs');
semilogy(areaSideLength_vec,nPLEstiError_ana2,'r--','LineWidth',1);

semilogy(areaSideLength_vec,nPLEstiError_simu3,'bo','MarkerFaceColor','b');
semilogy(areaSideLength_vec,nPLEstiError_ana3,'b-','LineWidth',1);

semilogy(areaSideLength_vec,nPLEstiError_simu4,'rs','MarkerFaceColor','r');
semilogy(areaSideLength_vec,nPLEstiError_ana4,'r-','LineWidth',1);


hold off;
xlabel('Square Area Side Length');
ylabel('$n_{\mathrm{PL}}$ Estimation Error','interpreter','latex');
legend('Grid $d=5$ (simu)','Grid $d=5$ (ana)', 'Random $\lambda=4\times 10^{-2}$ (simu)','Random $\lambda=4\times 10^{-2}$ (ana)',...
    'Grid $d=10$(simu)','Grid $d=10$(ana)', 'Random $\lambda=10^{-3}$ (simu)','Random $\lambda=10^{-3}$ (ana)','interpreter','latex');
grid on;
%xlim([50,500]);

figure;
semilogy(areaSideLength_vec,KdBEstiError_simu,'bo');
hold on;
semilogy(areaSideLength_vec,KdBEstiError_ana,'b--','LineWidth',1);
semilogy(areaSideLength_vec,KdBEstiError_simu2,'rs');
semilogy(areaSideLength_vec,KdBEstiError_ana2,'r--','LineWidth',1);

semilogy(areaSideLength_vec,KdBEstiError_simu3,'bo','MarkerFaceColor','b');
semilogy(areaSideLength_vec,KdBEstiError_ana3,'b-','LineWidth',1);

semilogy(areaSideLength_vec,KdBEstiError_simu4,'rs','MarkerFaceColor','r');
semilogy(areaSideLength_vec,KdBEstiError_ana4,'r-','LineWidth',1);


hold off;
xlabel('Square Area Side Length');
ylabel('$K_{\mathrm{dB}}$ Estimation Error','interpreter','latex');
legend('Grid $d=5$ (simu)','Grid $d=5$ (ana)', 'Random $\lambda=4\times 10^{-2}$ (simu)','Random $\lambda=4\times 10^{-2}$ (ana)',...
    'Grid $d=10$(simu)','Grid $d=10$(ana)', 'Random $\lambda=10^{-3}$ (simu)','Random $\lambda=10^{-3}$ (ana)','interpreter','latex');
grid on;
%xlim([50,300]);
