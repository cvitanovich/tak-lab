function [err,result] = errfun_specint5(V,IA_fullspec,IA_bandwidths)
%This function takes only ILDAlone mean arrays
%It attempts to combine bandwidth-delimited ILDAlone responses into a predicted
%broadband ILDAlone response while fitting coefficients

num_bands = size(IA_bandwidths,1);

%Specify indices for V
start_coeff_ind = 1; end_coeff_ind = start_coeff_ind - 1 + num_bands;

result = zeros(size(IA_fullspec));
for band_num = 1:num_bands
   result = result + ...
       (V(start_coeff_ind - 1 + band_num) *...
      (IA_bandwidths(band_num,:)));
end

index = find(result < 0);
result(index) = 0;

temp = corrcoef(IA_fullspec,result);
cc = temp(1,2);
err = sum(((result - IA_fullspec).^2))/sum((IA_fullspec).^2);
%err = 100 - (100*(cc.^2));

return