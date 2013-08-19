function err = ...
   errorfun_spectralintegration12(V,IA_fullspec,...
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
%   errorfun_spectralintegration10(V,IA_fullspec,...
%   ILD_piece1,...
%   ILD_piece2,...
%   ILD_piece3,...
%   ILD_piece4,...
%   ILD_piece5,...
%   ILD_piece6,...
%   ILD_piece7,...
%   ILD_piece8,...
%   IA_bestloc)
%V:               variables to fit
%IA_fullspec:     measured ILDAlone RS with broadband noise
%ILD_pieces:      ILD matrix from HRTF data (freqs assigned to rows, locs assigned to cols)
%                 divided into 8 1000-Hz wide bandwidths

%Assume that the maximal activity is broken up evenly among bandwidths
max_activity = max(IA_fullspec)/8;

%Assume that ILDwidths must be less than 20 dB
temp = V(17:24);
index_ILDwidth = find(temp > 10);
if(~isempty(index_ILDwidth))
   temp(index_ILDwidth) = 10;
   V(17:24) = temp;
end

%Assume that ILDopt_shift must be less than  +/- 10 dB
temp = V(9:16);
index_ILDopt_shift = find(temp > 10);
if(~isempty(index_ILDopt_shift))
   temp(index_ILDopt_shift) = 10;
   V(9:16) = temp;
end

%IA_piece1
ILD_opt_piece1 = repmat(ILD_piece1(:,IA_bestloc),1,size(ILD_piece1,2));
diff_ILD_piece1 = sqrt(mean((ILD_piece1 - ILD_opt_piece1 + V(9)).^2,1));
result_piece1 = max_activity * (1 - (diff_ILD_piece1/V(17)));

%IA_piece2
ILD_opt_piece2 = repmat(ILD_piece2(:,IA_bestloc),1,size(ILD_piece2,2));
diff_ILD_piece2 = sqrt(mean((ILD_piece2 - ILD_opt_piece2 + V(10)).^2,1));
result_piece2 = max_activity * (1 - (diff_ILD_piece2/V(18)));

%IA_piece3
ILD_opt_piece3 = repmat(ILD_piece3(:,IA_bestloc),1,size(ILD_piece3,2));
diff_ILD_piece3 = sqrt(mean((ILD_piece3 - ILD_opt_piece3 + V(11)).^2,1));
result_piece3 = max_activity * (1 - (diff_ILD_piece3/V(19)));

%IA_piece4
ILD_opt_piece4 = repmat(ILD_piece4(:,IA_bestloc),1,size(ILD_piece4,2));
diff_ILD_piece4 = sqrt(mean((ILD_piece4 - ILD_opt_piece4 + V(12)).^2,1));
result_piece4 = max_activity * (1 - (diff_ILD_piece4/V(20)));

%IA_piece5
ILD_opt_piece5 = repmat(ILD_piece5(:,IA_bestloc),1,size(ILD_piece5,2));
diff_ILD_piece5 = sqrt(mean((ILD_piece5 - ILD_opt_piece5 + V(13)).^2,1));
result_piece5 = max_activity * (1 - (diff_ILD_piece5/V(21)));

%IA_piece6
ILD_opt_piece6 = repmat(ILD_piece6(:,IA_bestloc),1,size(ILD_piece6,2));
diff_ILD_piece6 = sqrt(mean((ILD_piece6 - ILD_opt_piece6 + V(14)).^2,1));
result_piece6 = max_activity * (1 - (diff_ILD_piece6/V(22)));

%IA_piece7
ILD_opt_piece7 = repmat(ILD_piece7(:,IA_bestloc),1,size(ILD_piece7,2));
diff_ILD_piece7 = sqrt(mean((ILD_piece7 - ILD_opt_piece7 + V(15)).^2,1));
result_piece7 = max_activity * (1 - (diff_ILD_piece7/V(23)));

%IA_piece8
ILD_opt_piece8 = repmat(ILD_piece8(:,IA_bestloc),1,size(ILD_piece8,2));
diff_ILD_piece8 = sqrt(mean((ILD_piece8 - ILD_opt_piece8 + V(16)).^2,1));
result_piece8 = max_activity * (1 - (diff_ILD_piece8/V(24)));


result = (V(1)*result_piece1) + (V(2)*result_piece2)...
   + (V(3)*result_piece3) + (V(4)*result_piece4) + (V(5)*result_piece5)...
   + (V(6)*result_piece6) + (V(7)*result_piece7) + (V(8)*result_piece8);

temp = corrcoef(IA_fullspec,result);
cc = temp(1,2);
%err = sum(((result - IA_fullspec).^2))/sum((IA_fullspec).^2);
err = 100 - (100*(cc.^2));

return