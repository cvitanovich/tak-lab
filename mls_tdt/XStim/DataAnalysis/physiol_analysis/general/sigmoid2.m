function [y] = sigmoid2(x,height,beta,shiftx)
%Generate sigmoidal curve
%x 		: data vector
%height 	: maximal value of function
%beta 	: steepness of curve
%shiftx 	: x at half-height

if(nargin < 4) shiftx = 0; end

y = height./(1 + exp(-beta.*(x - shiftx)));

return