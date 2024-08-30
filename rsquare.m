function [R2 ci95] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model 
%
% [r2 ] = rsquare(y,f)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data y and model data f. 
% 
% INPUTS
%   y       : Actual data
%   f       : Model fit
%
%
% OUTPUT 
%   r2      : Coefficient of determination
% ci95:     : 95% confidenc interval



% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

R2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
n=length(y);
k=1;
SE = sqrt((4*R2*((1-R2).^2)*((n-k-1).^2))/(((n.^2)-1)*(n + 3)));
ci95=R2+[-1.96*SE 1.96*SE ];