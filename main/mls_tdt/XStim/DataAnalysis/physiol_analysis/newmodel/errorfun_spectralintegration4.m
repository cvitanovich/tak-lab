function [err,func] = errorfun_spectralintegration4(V,IA_fullspec,IA_bandwidths)
%This function takes only ILDAlone mean arrays
%It fits the following equation for data taken with 2 bandwidths:
%     result = (1/sum(V)) * (V(1) * IA_bandwidths(1,:)) + (V(2) * IA_bandwidths(2,:));

func = 'result = (1/sum(V)) * (V(1) * IA_bandwidths(1,:)) + (V(2) * IA_bandwidths(2,:)) + (V(3) * IA_bandwidths(3,:)) + (V(4) * IA_bandwidths(4,:));';

if (nargin == 0)
   err = -999;
   return
end

eval(func);

err = sum( (result - IA_fullspec).^2);

return