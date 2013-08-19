function [err,result,V] = lin_ilddistfun(V,IA_meansurf,dist_ild);

if(V(1) > 50) V(1) = 50; end
if(V(1) < 0) V(1) = 1; end
result = zeros(1,length(IA_meansurf));
result1 = max(IA_meansurf) * ...
   (1 - (abs(dist_ild(find(dist_ild<0)))/V(1)));
result2 = max(IA_meansurf) * ...
   (1 - (abs(dist_ild(find(dist_ild>0)))/V(1)));
result(find(dist_ild<0)) = result1;
result(find(dist_ild>0)) = result2;
result(find(dist_ild==0)) = max(IA_meansurf);
result(find(result < 0)) = 0;

err = sum( (IA_meansurf - result).^2)/sum(IA_meansurf.^2);

return