function [err,result,V] = gauss_ilddistfun(V,IA_meansurf,dist_ild);

if(V(1) > 50) V(1) = 50; end
if(V(1) < 0) V(1) = 1; end
result = max(IA_meansurf) * gaussian(dist_ild,0,V(1));

err = sum( (IA_meansurf - result).^2)/sum(IA_meansurf.^2);

return