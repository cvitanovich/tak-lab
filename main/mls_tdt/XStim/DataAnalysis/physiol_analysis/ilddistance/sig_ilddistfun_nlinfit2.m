function result = sig_ilddistfun_nlinfit2(beta,dist_ild,maxval);
%Order of beta:
%beta(1) : width or sigma of gaussian
%beta(2) : height of the sigmoidal part of the function
%beta(3) : beta for sigmoid2
%beta(4) : shift for sigmoid2

equation = normpdf(dist_ild,0,beta(1)) + sigmoid2(dist_ild,beta(2),beta(3),beta(4));
equation(find(equation < 0)) = 0;
equation = equation/max(equation);
result = maxval * equation;

return