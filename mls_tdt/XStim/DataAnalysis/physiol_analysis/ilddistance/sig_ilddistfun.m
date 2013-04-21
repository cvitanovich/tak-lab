function [err,result,V] = sig_ilddistfun(V,IA_meansurf,dist_ild);
%Order of V:
%V(1) : width or sigma of gaussian
%V(2) : height of the sigmoidal part of the function
%V(3) : beta for sigmoid2
%V(4) : shift for sigmoid2

if(V(1) > 50) V(1) = 50; end
if(V(1) < 4) V(1) = 4; end

equation = gaussian(dist_ild,0,V(1)) + sigmoid2(dist_ild,V(2),V(3),V(4));
equation(find(equation < 0)) = 0;
equation = equation/max(equation);
if(length(V) < 5)
   result = max(IA_meansurf) * equation;
else
   result = V(5) * equation;
end
err = sum( (IA_meansurf - result).^2)/sum(IA_meansurf.^2);

return