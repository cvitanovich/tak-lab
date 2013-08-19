function [err,result,V,freq_weight] = ...
   errfun_specint4(V,IA_fullspec,...
   if_surface,...
   if_freqaxis,...
   if_ildaxis,...
   trunc_ILD_matrix,...
   IA_bestloc_ind,...
   freq_bands,...
   freq_weight)
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
%Uses activity in ILD/freq along best ILDAlone location as frequency weights, instead of the
%activity peaks at each frequency

if(isempty(freq_bands))
   freq_bands(:,1:2) = [if_freqaxis if_freqaxis];
end

%specify indexing for V
start_ILDwidth = 1; end_ILDwidth = start_ILDwidth + size(freq_bands,1) - 1;
start_ILDshift = end_ILDwidth + 1; end_ILDshift = start_ILDshift + size(freq_bands,1) - 1;
lenV = length(V);

%Assume that ILDwidths must be less than a certain dB
max_allowable_ILDwidth = 25; %dB
temp = V(start_ILDwidth:end_ILDwidth);
index_ILDwidth = find(temp > max_allowable_ILDwidth);
if(~isempty(index_ILDwidth))
   temp(index_ILDwidth) = max_allowable_ILDwidth;
   V(start_ILDwidth:end_ILDwidth) = temp;
end

%Assume that ILDopt_shift must be less than  +/- some amount
max_allowable_ILDshift = 5; %dB
temp = V(start_ILDshift:end_ILDshift);
index_ILDshift = find(abs(temp) > max_allowable_ILDshift);
if(~isempty(index_ILDshift))
   temp(index_ILDshift) = sign(temp(index_ILDshift)) * max_allowable_ILDshift;
   V(start_ILDshift:end_ILDshift) = temp;
end

for band_num = 1:size(freq_bands,1)
   [y,minfreqindex] = min(abs(if_freqaxis - freq_bands(band_num,1)));
   [y,maxfreqindex] = min(abs(if_freqaxis - freq_bands(band_num,2)));
   
   ILD_opt = repmat(trunc_ILD_matrix(minfreqindex:maxfreqindex,IA_bestloc_ind),1,size(trunc_ILD_matrix,2));
   dist = sqrt(mean((trunc_ILD_matrix(minfreqindex:maxfreqindex,:) -...
      ILD_opt + V(start_ILDshift - 1 + band_num)).^2,1));
   result_pieces(band_num,:) = freq_weight(band_num) * ...
      (1 - ((dist/V(start_ILDwidth-1+band_num)).^V(lenV)));
end

result = sum(result_pieces,1);

index = find(result < 0);
result(index) = 0;

temp = corrcoef(IA_fullspec,result);
cc = temp(1,2);
%err = sum(((result - IA_fullspec).^2))/sum((IA_fullspec).^2);
err = 100 - (100*(cc.^2));

return
