function [err,result,V] = ...
   errfun_specint2(V,IA_fullspec,...
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

%specify indexing for V
start_ILDwidth = 1; end_ILDwidth = start_ILDwidth + size(freq_bands,1) - 1;
start_ILDshift = end_ILDwidth + 1; end_ILDshift = start_ILDshift + size(freq_bands,1) - 1;

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

%Limit frequencies according to freq_bands
min_freq = min(min(freq_bands));
max_freq = max(max(freq_bands));
minfreqind = min(find(if_freqaxis >= min_freq));
maxfreqind = max(find(if_freqaxis <= max_freq));

for freq_num = 1:length(if_freqaxis(minfreqind:maxfreqind))
   
   band_num = find(freq_bands(:,1) <= if_freqaxis(freq_num) & ...
      freq_bands(:,2) >= if_freqaxis(freq_num));
   band_num = max(band_num);
   
   ILD_opt = trunc_ILD_matrix(freq_num,IA_bestloc_ind);
   
   for loc_num = 1:size(trunc_ILD_matrix,2)
      [y,if_ildindex] = min(abs(if_ildaxis - trunc_ILD_matrix(freq_num,loc_num)));
      if_ildindex = max(if_ildindex);
      dist = sqrt(mean((trunc_ILD_matrix(freq_num,loc_num) -...
         ILD_opt + V(start_ILDshift - 1 + band_num)).^2,1));
      result_pieces(freq_num,loc_num) = if_surface(freq_num,if_ildindex) * (1 - (dist/V(start_ILDwidth-1+band_num)));
   end
end

result = sum(result_pieces,1);

index = find(result < 0);
result(index) = 0;

temp = corrcoef(IA_fullspec,result);
cc = temp(1,2);
%err = sum(((result - IA_fullspec).^2))/sum((IA_fullspec).^2);
err = 100 - (100*(cc.^2));

return
