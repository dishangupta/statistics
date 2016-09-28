alpha = 0.05;
zalphaby2 = norminv(1-alpha/2, 0 , 1); 

a = 0;
b = 1;
n = 100;

theta = [-3:0.01:3];
CovCn = zeros(size(theta));
CovDn = zeros(size(theta));

for i = 1:length(theta)
    meanCn = (theta(i) - a)/(b*sqrt(1+n*b^2));
    varCn = n*b^2/(1+n*b^2);
    CovCn(i) = normcdf(zalphaby2, meanCn, varCn) - normcdf(-zalphaby2, meanCn, varCn);
    
    meanDn = 0;
    varDn = 1;
    CovDn(i) = normcdf(zalphaby2, meanDn, varDn) - normcdf(-zalphaby2, meanDn, varDn); 
end

plot(theta, CovCn, theta, CovDn);