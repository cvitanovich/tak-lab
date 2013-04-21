function err = ...
   errorfun_spectralintegration7(V,IA_fullspec,IA_piece1,IA_piece2,ILD_piece1,ILD_piece2,IA_bestloc)
%This function takes only ILDAlone mean arrays
%errorfun_spectralintegration7(V,IA_fullspec,IA_piece1,IA_piece2,ILD_piece1,ILD_piece2,IA_bestloc)
%V:               variables to fit
%IA_fullspec:     measured ILDAlone RS with broadband noise
%IA_pieces:       bandwidth-limited measured ILDAlone RS's
%ILD_pieces:      ILD matrix from HRTF data (freqs assigned to rows, locs assigned to cols)

%IA_piece1
maxactivity_piece1 = max(IA_piece1);
ILD_opt_piece1 = repmat(ILD_piece1(:,IA_bestloc),1,size(ILD_piece1,2));
diff_ILD_piece1 = sqrt(mean((ILD_piece1 - ILD_opt_piece1).^2,1));
result_piece1 = maxactivity_piece1 * (1 - (diff_ILD_piece1/V(1)));

%IA_piece2
maxactivity_piece2 = max(IA_piece2);
ILD_opt_piece2 = repmat(ILD_piece2(:,IA_bestloc),1,size(ILD_piece2,2));
diff_ILD_piece2 = sqrt(mean((ILD_piece2 - ILD_opt_piece2).^2,1));
result_piece2 = maxactivity_piece2 * (1 - (diff_ILD_piece2/V(2)));

result = result_piece1 + result_piece2;
err = mean( (result - IA_fullspec).^2);

return