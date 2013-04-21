function [err,result,V] = sig_ilddistfun2(V,IA_meansurf,dist_ild);
%Order of V:
%V(1) : width or sigma of gaussian
%V(2) : height of the sigmoidal part of the function
%V(3) : beta for sigmoid2
%V(4) : shift for sigmoid2

if(V(1) > 50) V(1) = 50; end
if(V(1) < 0) V(1) = 1; end

equation = sigmoid2(dist_ild,V(2),V(3),V(4)) + sigmoid2(dist_ild,V(6),V(7),V(8));
equation(find(equation < 0)) = 0;
equation = equation/max(equation);
result = V(5) * equation;
err = sum( (IA_meansurf - result).^2)/sum(IA_meansurf.^2);

return