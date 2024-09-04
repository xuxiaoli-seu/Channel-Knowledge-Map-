%This file is created by Xu Xiaoli on 07/10/2023
%It estimates the shadowing parameters for based on the LS methods

%alpha: shadowing variance
%beta: correlation distance
%sigma: multipath fading variance

%YQ_shaFad: the channel gain substract the path loss, equals to
%shadowing+fading

function [alpha, beta,sigma_sq]=EstShadowPara_LS(YQ_shaFad,MeasureLoc)


k=length(YQ_shaFad); % the number of measurements
distMatrix=sqrt((repmat(MeasureLoc(1,:),k,1)-repmat(MeasureLoc(1,:)',1,k)).^2....
    +(repmat(MeasureLoc(2,:),k,1)-repmat(MeasureLoc(2,:)',1,k)).^2);

Chi_esti=(YQ_shaFad'*YQ_shaFad)/k; %the total variance, equals alpha+sigma^2

distAll=unique(distMatrix);

distThr=find(distAll>100,1,'first'); %distance beyond threshold is not considered 
if isempty(distThr)
    distThr=length(distAll);
end

distVec=distAll(2:distThr); %the niminimum distance is 0, hence we start from index 2
%Count the number of pairs with the distance
Yqi=cell(1,length(distVec));
Yqj=cell(1,length(distVec));

numDistance=zeros(1,length(distVec));
rQ_vec=zeros(1,length(distVec));
for j=1:length(distVec)
    for i=1:k-1
        indx=find(abs(distMatrix(i,i:k)-distVec(j))<0.5);
        Yqi{j}=[Yqi{j};repmat(YQ_shaFad(i),length(indx),1)];
        Yqj{j}=[Yqj{j};YQ_shaFad(indx+i-1)];
    end
    numDistance(j)=length(Yqi{j});
    measure1=Yqi{j};
    measure2=Yqj{j};
    rQ_vec(j)=measure1'*measure2/numDistance(j);
    %If negative correlation is observed, it implies that almost no
    %correlation beyond this point, and the rest of data is neglected
    dist_thre=j;
    if rQ_vec(j)<0
        rQ_vec(j)=0.0001; 
        break;
    end 
end

MLQ=[ones(dist_thre,1),-distVec(1:dist_thre)];
b=log(rQ_vec(1:dist_thre)');
WLQ=diag(numDistance(1:dist_thre)); %change the weight for correlation distance estimation

para=(MLQ/(MLQ'*WLQ*MLQ)')'*WLQ*b;
alpha=exp(para(1));
beta=1/para(2);

if (alpha<Chi_esti)&&(beta>0)
    sigma_sq=Chi_esti-alpha;
elseif (alpha>Chi_esti &&(beta>0))
    alpha=Chi_esti;sigma_sq=0;
    beta=(numDistance(1:dist_thre)*(distVec(1:dist_thre)./(log(alpha)-b)))/sum(numDistance(1:dist_thre));
else
    alpha=0;beta=0;sigma_sq=Chi_esti;
end
%sigma_sq=Chi_esti-alpha;