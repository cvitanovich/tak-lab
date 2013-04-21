function [err,result] = errfun_specint6(V,IA_fullspec,IA_bandwidths)
%This function takes only ILDAlone mean arrays
%It attempts to combine bandwidth-delimited ILDAlone responses into a predicted
%broadband ILDAlone response by using a threshhold

num_bands = size(IA_bandwidths,1);
maxIAfullspec = max(IA_fullspec);

result = zeros(size(IA_fullspec));
for band_num = 1:num_bands
   thresh_ind = find(IA_bandwidths(band_num,:) < (V(2*band_num)*maxIAfullspec));
   IA_bandwidths(band_num,thresh_ind) = 0;
   result = result + ...
       (V(band_num)) .*...
      (IA_bandwidths(band_num,:));
end

index = find(result < 0);
result(index) = 0;

temp = corrcoef(IA_fullspec,result);
cc = temp(1,2);
%err = sum(((result - IA_fullspec).^2))/sum((IA_fullspec).^2);
err = 100 - (100*(cc.^2));

return