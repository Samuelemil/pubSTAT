function p=aucCompareIndependent(auc1,nPositive1,nNegative1,auc2,nPositive2,nNegative2)
% Test if two ROC curves are significant difference in two Independent
% datasets (two sided, in current implementation)
% 
%  p=aucCompareIndependent(A1,nPositive1,nNegative1,A2,nPositive2,nNegative2)
% auc1:           AUC in dataset 1
% nPositive1:   n positive cases in dataset 1
% nNegative1:   n negative cases in dataset 1
% auc2:           AUC in dataset 3
% nPositive2:   n positive cases in dataset 2
% nNegative2:   n negative cases in dataset 2
%
% p:            p valuse two sided 
% Hanley JA, McNeil BJ (1982) The meaning and use of the area under a receiver operating characteristic (ROC) curve. Radiology 143:29-36.
%  Author: Samuel E Schmidt sschmidt@hst.aau.dk



A=auc1;
Pxxy=A/(2-A);
Pxyy=(2*A^2)/(1+A);
se_a=((A*(1-A)+(nPositive1-1)*(Pxxy-A^2)+(nNegative1-1)*(Pxyy-A^2))/(nPositive1*nNegative1))^(1/2);


A=auc2;
Pxxy=A/(2-A);
Pxyy=(2*A^2)/(1+A);
se_b=((A*(1-A)+(nPositive2-1)*(Pxxy-A^2)+(nNegative2-1)*(Pxyy-A^2))/(nPositive2*nNegative2))^(1/2);
se_sum=sqrt(se_b^2+se_a^2);

A_diff=auc1-auc2;

z=abs(A_diff/se_sum);

 p=2*(1-normcdf(z,0,1));
