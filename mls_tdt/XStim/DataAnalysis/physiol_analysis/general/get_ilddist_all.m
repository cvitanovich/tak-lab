function [dist_ild] = get_ilddist_all(IA_HRTFmat,IA_bestloc_ind);
%Inputs:
%	IA_HRTFmat for the corresponding locations
%	ILDAlone_bestlocation_index
%  Output the matrix quantifying distance from optimal ILDAlone ILD spectrum

IA_HRTFmat_bestloc = repmat(IA_HRTFmat(:,IA_bestloc_ind),1,size(IA_HRTFmat,2));
dist_ild = sqrt((IA_HRTFmat - IA_HRTFmat_bestloc).^2);
return