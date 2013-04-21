function [err,Tonal_IA_meansurf] = ilddist_weights(V,ILDf_mat,...
   ILD_matrix,...
   HRTF_freqaxis,...
   if_meansurf,...
   if_freqaxis,...
   if_ildaxis,...
   IA_bestloc_ind,...
   IA_meansurf);
%Function that implements Line Integral method WITH ILD_distance correction
%V:					A Vector of ILDwidths for minimization of the function
%ILDf_mat_all:		A matrix of ILDxfreq matrices, 1 for each location
%ILD_matrix:		ILD matrix (r=freqs,c=locs)
%HRTF_freqaxis:	frequencies for ILD matrix
%if_meansurf:		ILD/freq matrix (r=ILDs,c=freqs)
%if_freqaxis:		frequencies for the ILD/freq matrix
%if_ildaxis:		ILDs for the ILD/freq matrix
%IA_bestloc_ind:	index for the location of maximal activity in the ILDAlone RS

maxILD = max(if_ildaxis); minILD = min(if_ildaxis);

for freq_num = 1:length(if_freqaxis)
   [y,hrtf_freqind] = min(abs(HRTF_freqaxis - if_freqaxis(freq_num)));
   dist_ild = abs(ILD_matrix(hrtf_freqind,:) - ...
      ILD_matrix(hrtf_freqind,IA_bestloc_ind));
   tempmat(freq_num,:) = ILDf_mat(freq_num,:) .* ...
      (1 - ((dist_ild)/V(freq_num)));
end

Tonal_IA_meansurf = sum(tempmat,1);
neg_index = find(Tonal_IA_meansurf < 0);
Tonal_IA_meansurf(neg_index) = 0;

temp = corrcoef(IA_meansurf,Tonal_IA_meansurf);
vari_exp = temp(1,2)^2;
err = 100*(1 - vari_exp);

%disp('Evaluated function...')
         
   