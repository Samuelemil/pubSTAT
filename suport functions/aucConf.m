function conf=aucConf(auc,nA,nN, oneOrTwoSided)
% Estimates the 95% confidence interval as either one sides (the lower interval) or two sided. 
% aucConf(auc,nA,nN, oneSided)
% auc: esitmated auc value
% nA,number of  Diseased subjects
% nN: number of Healthy subjects 
% oneOrTwoSided: One sided (=1) or two sided test (=2)
% See Hanley, J.A.,The meaning and use of the area under a receiver operating characteristic (ROC) curve, 1982
%  Author: Samuel E Schmidt sschmidt@hst.aau.dk

% estimate parameters
Q1=auc/(2-auc);
Q2=2*(auc^2)/(1+auc);

% estimate standart error
se=((auc*(1-auc)+(nA-1)*(Q1-auc^2)+(nN-1)*(Q2-auc^2))/(nA*nN))^(1/2);


% estiamte confidence intervals
if oneOrTwoSided==1
z=1.645;
conf=[auc-z*se];
elseif oneOrTwoSided==2
   z=1.96; 

   conf= [auc-z*se auc+z*se];
end

