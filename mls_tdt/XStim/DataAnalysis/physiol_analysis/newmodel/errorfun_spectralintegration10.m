function err = ...
   errorfun_spectralintegration10(V,IA_fullspec,...
   ILD_piece1,...
   ILD_piece2,...
   ILD_piece3,...
   ILD_piece4,...
   ILD_piece5,...
   ILD_piece6,...
   ILD_piece7,...
   ILD_piece8,...
   IA_bestloc)
%This function takes only ILDAlone mean arrays
%errorfun_spectralintegration7(V,IA_fullspec,IA_piece1,IA_piece2,ILD_piece1,ILD_piece2,IA_bestloc)
%V:               variables to fit
%IA_fullspec:     measured ILDAlone RS with broadband noise
%IA_pieces:       bandwidth-limited measured ILDAlone RS's
%ILD_pieces:      ILD matrix from HRTF data (freqs assigned to rows, locs assigned to cols)

max_activity = max(IA_fullspec)/8;

%IA_piece1
ILD_opt_piece1 = repmat(ILD_piece1(:,IA_bestloc),1,size(ILD_piece1,2));
diff_ILD_piece1 = sqrt(mean((ILD_piece1 - ILD_opt_piece1).^2,1));
result_piece1 = max_activity * (1 - (diff_ILD_piece1/V(9)));

%IA_piece2
ILD_opt_piece2 = repmat(ILD_piece2(:,IA_bestloc),1,size(ILD_piece2,2));
diff_ILD_piece2 = sqrt(mean((ILD_piece2 - ILD_opt_piece2).^2,1));
result_piece2 = max_activity * (1 - (diff_ILD_piece2/V(9)));

%IA_piece3
ILD_opt_piece3 = repmat(ILD_piece3(:,IA_bestloc),1,size(ILD_piece3,2));
diff_ILD_piece3 = sqrt(mean((ILD_piece3 - ILD_opt_piece3).^2,1));
result_piece3 = max_activity * (1 - (diff_ILD_piece3/V(9)));

%IA_piece4
ILD_opt_piece4 = repmat(ILD_piece4(:,IA_bestloc),1,size(ILD_piece4,2));
diff_ILD_piece4 = sqrt(mean((ILD_piece4 - ILD_opt_piece4).^2,1));
result_piece4 = max_activity * (1 - (diff_ILD_piece4/V(9)));

%IA_piece5
ILD_opt_piece5 = repmat(ILD_piece5(:,IA_bestloc),1,size(ILD_piece5,2));
diff_ILD_piece5 = sqrt(mean((ILD_piece5 - ILD_opt_piece5).^2,1));
result_piece5 = max_activity * (1 - (diff_ILD_piece5/V(10)));

%IA_piece6
ILD_opt_piece6 = repmat(ILD_piece6(:,IA_bestloc),1,size(ILD_piece6,2));
diff_ILD_piece6 = sqrt(mean((ILD_piece6 - ILD_opt_piece6).^2,1));
result_piece6 = max_activity * (1 - (diff_ILD_piece6/V(10)));

%IA_piece7
ILD_opt_piece7 = repmat(ILD_piece7(:,IA_bestloc),1,size(ILD_piece7,2));
diff_ILD_piece7 = sqrt(mean((ILD_piece7 - ILD_opt_piece7).^2,1));
result_piece7 = max_activity * (1 - (diff_ILD_piece7/V(10)));

%IA_piece8
ILD_opt_piece8 = repmat(ILD_piece8(:,IA_bestloc),1,size(ILD_piece8,2));
diff_ILD_piece8 = sqrt(mean((ILD_piece8 - ILD_opt_piece8).^2,1));
result_piece8 = max_activity * (1 - (diff_ILD_piece8/V(10)));


result = (V(1)*result_piece1) + (V(2)*result_piece2)...
   + (V(3)*result_piece3) + (V(4)*result_piece4) + (V(5)*result_piece5)...
   + (V(6)*result_piece6) + (V(7)*result_piece7) + (V(8)*result_piece8);

%err = sum(((result - IA_fullspec).^2))/sum((IA_fullspec).^2);
temp = corrcoef(IA_fullspec,result);
cc = temp(1,2);
err = 100 - (100*(cc.^2));

return