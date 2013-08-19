function [Tonal_IA_meansurf] = if2space2(ILD_matrix,...
   HRTF_freqaxis,...
   if_meansurf,...
   if_freqaxis,...
   if_ildaxis,...
   IA_bestloc_ind);
%Function that implements Line Integral method WITH ILD_distance correction
%ILD_matrix:		ILD matrix (r=freqs,c=locs)
%HRTF_freqaxis:	frequencies for ILD matrix
%if_meansurf:		ILD/freq matrix (r=ILDs,c=freqs)
%if_freqaxis:		frequencies for the ILD/freq matrix
%if_ildaxis:		ILDs for the ILD/freq matrix
%IA_bestloc_ind:	index for the location of maximal activity in the ILDAlone RS

%Delimit frequencies
minfreq = min(if_freqaxis); maxfreq = max(if_freqaxis);
[y,minfreqind] = min(abs(HRTF_freqaxis - minfreq));
[y,maxfreqind] = min(abs(HRTF_freqaxis - maxfreq));

ILD_matrix = ILD_matrix(minfreqind:maxfreqind,:);
HRTF_freqaxis = HRTF_freqaxis(minfreqind:maxfreqind);

%Get ILD limits
minILD = min(if_ildaxis); maxILD = max(if_ildaxis);

for loc_num = 1:size(ILD_matrix,2)
   %Transform spectrum into a matrix of 0's and 1's
   ILDf_mat = zeros(size(if_meansurf,1),size(if_meansurf,2));
   for freq_num = 1:length(if_freqaxis)
      [y,freqind] = min(abs(HRTF_freqaxis - if_freqaxis(freq_num)));
      [y,ildind]  = min(abs(if_ildaxis - ILD_matrix(freqind,loc_num)));
      %Make sure the ILD value from ILD_matrix is in the range of the ILDfreq ildaxis
      if(ILD_matrix(freqind,loc_num) <= maxILD & ILD_matrix(freqind,loc_num) >= minILD)
         ILDf_mat(ildind,freq_num) = 1;
      end
   end
   temp = if_meansurf .* ILDf_mat;
   Tonal_IA_meansurf(loc_num) = sum(temp(:));
   ind_neg = find(Tonal_IA_meansurf < 0);
   Tonal_IA_meansurf(ind_neg) = 0;
   %disp(['Finished location ' num2str(loc_num)])
end