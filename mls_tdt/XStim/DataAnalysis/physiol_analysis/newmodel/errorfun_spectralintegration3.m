function [err,func] = errorfun_spectralintegration3(V,IA_fullspec,IA_bandwidths)
%This function takes only ILDAlone mean arrays
%It fits the following equation for data taken with 2 bandwidths:
%     result = (1/sum(V)) * (V(1) * IA_bandwidths(1,:)) + (V(2) * IA_bandwidths(2,:));

func = 'result = ((V(1) * IA_bandwidths(1,:)) + (V(2) * IA_bandwidths(2,:)));';

if (nargin == 0)
   err = -999;
   return
end

eval(func);

err = mean( (result - IA_fullspec).^2);

return