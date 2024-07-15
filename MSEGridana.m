function MSE=MSEGridana(d,k_vec,alpha,beta,sigma)
    var=d/beta;
    func=exp(-var)*(0.2854-0.0725*var+0.0108*var^2);
    para=pi*beta^2/(2*d^2)*(1-(1+sqrt(2)*var)*exp(-sqrt(2)*var))-2*func;
    MSE=alpha+sigma^2-(k_vec*alpha^2)./(k_vec*alpha+sigma^2)*para;   
end