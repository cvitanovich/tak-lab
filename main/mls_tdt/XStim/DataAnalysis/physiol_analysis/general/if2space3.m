function [ILDf_mat_all] = if2space3(ILD_matrix,...
   HRTF_freqaxis,...
   if_meansurf,...
   if_freqaxis,...
   if_ildaxis);
%Function that implements Line Integral method WITH ILD_distance correction
%ILD_matrix:		ILD matrix (r=freqs,c=locs)
%HRTF_freqaxis:	frequencies for ILD matrix
%if_meansurf:		ILD/freq matrix (r=ILDs,c=freqs)
%if_freqaxis:		frequencies for the ILD/freq matrix
%if_ildaxis:		ILDs for the ILD/freq matrix

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
   ILDf_mat_all(:,:,loc_num) = temp;
end