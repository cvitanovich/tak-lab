function [ild_bands] = get_ild_bands(IA_HRTFmat,...
   IA_HRTFmat_freqs,...
   freq_bands,...
   sign_flag);
%Inputs:
%	IA_HRTFmat for the corresponding locations
%  IA_HRTFmat_freqs frequencies for the matrix
%	ILDAlone_bestlocation_index
%	the freq limits, [min max], (n x 2)
%  sign_flag: 0, all values positive; 1, values signed
%Output the axis quantifying distance from optimal ILDAlone ILD spectrum

for flim_num = 1:size(freq_bands,1)
   [y,minfreqind] = min(abs(IA_HRTFmat_freqs - freq_bands(flim_num,1)));
   [y,maxfreqind] = min(abs(IA_HRTFmat_freqs - freq_bands(flim_num,2)));
   if(sign_flag)
      ild_bands(flim_num,:) = mean(IA_HRTFmat(minfreqind:maxfreqind,:),1);
   else
      ild_bands(flim_num,:) = sqrt(mean(IA_HRTFmat(minfreqind:maxfreqind,:).^2,1));
   end   
end

return