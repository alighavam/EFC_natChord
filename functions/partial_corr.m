function [rho,p,resX,resY] = partial_corr(X,Y,Z)
% regress out Z from X and Y and calculate correlation of residuals.

beta1 = svd_linregress(Y,Z);
resY = Y - Z*beta1;

beta2 = svd_linregress(X,Z);
resX = X - Z*beta2;

rho = corr(resY,resX);

mdl = fitlm(resY,resX);
p = mdl.Coefficients.pValue(2);




