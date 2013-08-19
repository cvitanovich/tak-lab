function [err,result,V,freq_weight] = ...
   errfun_specint8(V,IA_fullspec,...
   if_surface,...
   if_freqaxis,...
   if_ildaxis,...
   trunc_ILD_matrix,...
   IA_bestloc_ind,...
   freq_bands)
%error function for spectral integration
%Each location's contribution is weighted by
%		1) in frequency by that frequency's Tonal Response
%		2) in ILD by the distance from from the optimal ILD spectrum
%if_surface:	ILD/freq surface (tonal or bandpass) *after* interpolation
%if_freqaxis:	Frequency axis of if_surface
%if_ildaxis:	ILD axis of if_surface
%trunc_ILD_matrix:	ILD_matrix containing only those freqs in if_surface
%IA_bestloc_ind:	index of the best ILDAlone location in the ILD_matrix
%freq_bands:	a Nx2 matrix containing N [minfreq maxfreq] rows of the frequency bands

if(isempty(freq_bands))
   freq_bands(:,1:2) = [if_freqaxis if_freqaxis];
end

lenV = length(V);

tempind = find(V < 0);
V(tempind) = 1;
tempind = find(V > 30);
V(tempind) = 30;

maxval = max(IA_fullspec);

for band_num = 1:size(freq_bands,1)
   [y,minfreqindex] = min(abs(if_freqaxis - freq_bands(band_num,1)));
   [y,maxfreqindex] = min(abs(if_freqaxis - freq_bands(band_num,2)));
   
   freq_weight(band_num) = mean(max(if_surface(minfreqindex:maxfreqindex,:),[],1));
   
   ILD_opt = repmat(trunc_ILD_matrix(minfreqindex:maxfreqindex,IA_bestloc_ind),1,size(trunc_ILD_matrix,2));
   dist = sqrt(mean((trunc_ILD_matrix(minfreqindex:maxfreqindex,:) -...
      ILD_opt).^2,1));
   P = (1 - ((dist/V(band_num))));
   if(P>=0)
      result_pieces(band_num,:) = (freq_weight(band_num) * maxval ...
         * (P.^V(lenV)))./(V(lenV-1).^V(lenV) + P.^V(lenV));
   else
      result_pieces(band_num,:) = 0 * P;
   end
end

result = sum(result_pieces,1);

index = find(result < 0);
result(index) = 0;

temp = corrcoef(IA_fullspec,result);
cc = temp(1,2);
err = 100 - (100*(cc.^2));
%err = sum(((result - IA_fullspec).^2))/sum((IA_fullspec).^2);

return
