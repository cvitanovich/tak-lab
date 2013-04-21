function [dist_ild] = get_ilddist(IA_HRTFmat,...
   IA_HRTFmat_freqs,...
   IA_bestloc_ind,...
   freq_bands,...
   sign_flag);
%Inputs:
%	IA_HRTFmat for the corresponding locations
%  IA_HRTFmat_freqs frequencies for the matrix
%	ILDAlone_bestlocation_index
%	the freq limits, [min max], (n x 2)
%  sign_flag: 0, all values positive; 1, values signed
%Output the axis quantifying distance from optimal ILDAlone ILD spectrum

disp('Getting ILD distance!')
for flim_num = 1:size(freq_bands,1)
   [y,minfreqind] = min(abs(IA_HRTFmat_freqs - freq_bands(flim_num,1)));
   [y,maxfreqind] = min(abs(IA_HRTFmat_freqs - freq_bands(flim_num,2)));
   IA_HRTFmat_bestloc = repmat(IA_HRTFmat(minfreqind:maxfreqind,IA_bestloc_ind),1,size(IA_HRTFmat,2));
   if(sign_flag)
      dist_ild(flim_num,:) = sqrt(mean((IA_HRTFmat(minfreqind:maxfreqind,:) - ...
         IA_HRTFmat_bestloc).^2,1)) .* ...
      sign(mean(IA_HRTFmat(minfreqind:maxfreqind,:) - IA_HRTFmat_bestloc,1));
   else
      dist_ild(flim_num,:) = sqrt(mean((IA_HRTFmat(minfreqind:maxfreqind,:) - ...
         IA_HRTFmat_bestloc).^2,1));
   end   
end

return