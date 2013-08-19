function [dist_ild] = get_ilddist2(HRTFmat,HRTFmat_freqs,IA_HRTFbestloc, freq_bands,sign_flag);

%function [dist_ild] = get_ilddist2(HRTFmat,HRTFmat_freqs,IA_HRTFbestloc,f_bands,sign_flag);
%
%Inputs:
%	HRTFmat (ILD in dB for each location as nfreqs x nLocs)
%  HRTFmat_freqs frequencies for the matrix
%	DIFFERENCE FROM VER 1.0:  ILD from ILDAlone_bestlocation  (not just the index, but actually the ILD)
%	the freqband limits in Hz, [minf maxf] as (n x 2)
%  sign_flag: 0, all values positive; 1, values signed
%Output the axis quantifying distance from optimal ILDAlone ILD spectrum

%disp('Getting ILD distance!')
IA_HRTFmatbestloc = repmat(IA_HRTFbestloc(:),1,size(HRTFmat,2));

for flim_num = 1:size(freq_bands,1)
   [y,minfreqind] = min(abs(HRTFmat_freqs - freq_bands(flim_num,1)));
   [y,maxfreqind] = min(abs(HRTFmat_freqs - freq_bands(flim_num,2)));
   IA_HRTFmat_bestloc = repmat(IA_HRTFbestloc(minfreqind:maxfreqind),1,size(IA_HRTFmatbestloc,2));
   if(sign_flag)
      dist_ild(flim_num,:) = sqrt(mean((HRTFmat(minfreqind:maxfreqind,:) - ...
         IA_HRTFmat_bestloc).^2,1)) .* sign(mean(HRTFmat(minfreqind:maxfreqind,:) - IA_HRTFmat_bestloc,1));
   else
      dist_ild(flim_num,:) = sqrt(mean((HRTFmat(minfreqind:maxfreqind,:) - IA_HRTFmat_bestloc).^2,1));
   end   
end

return
