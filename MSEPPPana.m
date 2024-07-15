function MSE=MSEPPPana(density,k_vec,alpha,beta,sigma)
   
    coef=pi*density*beta^2;
    para= (1-exp(1./coef)./(beta*sqrt(density)).*(1-erf(1./sqrt(coef))));
    
    MSE=alpha+sigma^2-k_vec.*alpha^2./(k_vec*alpha+sigma^2)*para;