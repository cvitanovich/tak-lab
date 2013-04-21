function err = ...
   errorfun_spectralintegration13(V,IA_fullspec,...
   ILD_piece1,...
   ILD_piece2,...
   IA_bestloc)
%This function takes only ILDAlone mean arrays
%   errorfun_spectralintegration10(V,IA_fullspec,...
%   ILD_piece1,...
%   ILD_piece2,...
%   IA_bestloc)
%V:               variables to fit
%IA_fullspec:     measured ILDAlone RS with broadband noise
%ILD_pieces:      ILD matrix from HRTF data (freqs assigned to rows, locs assigned to cols)
%                 divided into 8 1000-Hz wide bandwidths

%Assume that the maximal activity is broken up evenly among bandwidths
max_activity = max(IA_fullspec)/8;

%Assume that weights must be >= 0
temp = V(1:2);
index_weights = find(temp < 0);
if(~isempty(index_weights))
   temp(index_weights) = 0;
   V(1:2) = temp;
end

%Assume that ILDwidths must be less than 20 dB
temp = V(5:6);
index_ILDwidth = find(temp > 10);
if(~isempty(index_ILDwidth))
   temp(index_ILDwidth) = 10;
   V(5:6) = temp;
end

%Assume that ILDopt_shift must be less than  +/- 10 dB
temp = V(3:4);
index_ILDopt_shift = find(temp > 10);
if(~isempty(index_ILDopt_shift))
   temp(index_ILDopt_shift) = 10;
   V(3:4) = temp;
end

%IA_piece1
ILD_opt_piece1 = repmat(ILD_piece1(:,IA_bestloc),1,size(ILD_piece1,2));
diff_ILD_piece1 = sqrt(mean((ILD_piece1 - ILD_opt_piece1 + V(3)).^2,1));
result_piece1 = max_activity * (1 - (diff_ILD_piece1/V(5)));

%IA_piece2
ILD_opt_piece2 = repmat(ILD_piece2(:,IA_bestloc),1,size(ILD_piece2,2));
diff_ILD_piece2 = sqrt(mean((ILD_piece2 - ILD_opt_piece2 + V(4)).^2,1));
result_piece2 = max_activity * (1 - (diff_ILD_piece2/V(6)));

result = (V(1)*result_piece1) + (V(2)*result_piece2);
index_neg = find(result < 0);
result(index_neg) = 0;


temp = corrcoef(IA_fullspec,result);
cc = temp(1,2);
%err = sum(((result - IA_fullspec).^2))/sum((IA_fullspec).^2);
err = 100 - (100*(cc.^2));

return