function err = ...
   errorfun_spectralintegration14(V,IA_fullspec,...
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

%specify indexing for V
start_weights = 1; end_weights = 8;
start_ILDwidth = 9; end_ILDwidth = 16;
start_ILDpeakshift1 = 17; end_ILDpeakshift1 = start_ILDpeakshift1 + size(ILD_piece1,1)-1;
start_ILDpeakshift2 = end_ILDpeakshift1 + 1; end_ILDpeakshift2 = start_ILDpeakshift2 + size(ILD_piece2,1)-1;
start_ILDpeakshift3 = end_ILDpeakshift2 + 1; end_ILDpeakshift3 = start_ILDpeakshift3 + size(ILD_piece3,1)-1;
start_ILDpeakshift4 = end_ILDpeakshift3 + 1; end_ILDpeakshift4 = start_ILDpeakshift4 + size(ILD_piece4,1)-1;
start_ILDpeakshift5 = end_ILDpeakshift4 + 1; end_ILDpeakshift5 = start_ILDpeakshift5 + size(ILD_piece5,1)-1;
start_ILDpeakshift6 = end_ILDpeakshift5 + 1; end_ILDpeakshift6 = start_ILDpeakshift6 + size(ILD_piece6,1)-1;
start_ILDpeakshift7 = end_ILDpeakshift6 + 1; end_ILDpeakshift7 = start_ILDpeakshift7 + size(ILD_piece7,1)-1;
start_ILDpeakshift8 = end_ILDpeakshift7 + 1; end_ILDpeakshift8 = start_ILDpeakshift8 + size(ILD_piece8,1)-1;

%Assume that the maximal activity is broken up evenly among bandwidths
max_activity = max(IA_fullspec)/8;

%Assume that ILDwidths must be less than 20 dB
temp = V(start_ILDwidth:end_ILDwidth);
index_ILDwidth = find(temp > 10);
if(~isempty(index_ILDwidth))
   temp(index_ILDwidth) = 10;
   V(start_ILDwidth:end_ILDwidth) = temp;
end

%Assume that ILDopt_shift must be less than  +/- 5 dB
num_pieces = 8;
max_ILDopt_shift = 5;
for n_piece = 1:num_pieces
   eval(['start_ind = start_ILDpeakshift' num2str(n_piece) ';']);
   eval(['end_ind = end_ILDpeakshift' num2str(n_piece) ';']);
   temp = V(start_ind:end_ind);
   index_ILDopt_shift = find(abs(temp) > max_ILDopt_shift);
   if(~isempty(index_ILDopt_shift))
      temp(index_ILDopt_shift) = max_ILDopt_shift * sign(temp(index_ILDopt_shift));
      V(start_ind:end_ind) = temp;
   end
   clear temp start_ind end_ind
end

%IA_piece1
V_mat = repmat(V(start_ILDpeakshift1:end_ILDpeakshift1)',1,size(ILD_piece1,2));
ILD_opt_piece1 = repmat(ILD_piece1(:,IA_bestloc),1,size(ILD_piece1,2));
diff_ILD_piece1 = sqrt(mean((ILD_piece1 - ILD_opt_piece1 + V_mat).^2,1));
result_piece1 = max_activity * (1 - (diff_ILD_piece1/V(start_ILDwidth)));

%IA_piece2
V_mat = repmat(V(start_ILDpeakshift2:end_ILDpeakshift2)',1,size(ILD_piece2,2));
ILD_opt_piece2 = repmat(ILD_piece2(:,IA_bestloc),1,size(ILD_piece2,2));
diff_ILD_piece2 = sqrt(mean((ILD_piece2 - ILD_opt_piece2 + V_mat).^2,1));
result_piece2 = max_activity * (1 - (diff_ILD_piece2/V(start_ILDwidth+1)));

%IA_piece3
V_mat = repmat(V(start_ILDpeakshift3:end_ILDpeakshift3)',1,size(ILD_piece2,2));
ILD_opt_piece3 = repmat(ILD_piece3(:,IA_bestloc),1,size(ILD_piece3,2));
diff_ILD_piece3 = sqrt(mean((ILD_piece3 - ILD_opt_piece3 + V_mat).^2,1));
result_piece3 = max_activity * (1 - (diff_ILD_piece3/V(start_ILDwidth+2)));

%IA_piece4
V_mat = repmat(V(start_ILDpeakshift4:end_ILDpeakshift4)',1,size(ILD_piece2,2));
ILD_opt_piece4 = repmat(ILD_piece4(:,IA_bestloc),1,size(ILD_piece4,2));
diff_ILD_piece4 = sqrt(mean((ILD_piece4 - ILD_opt_piece4 + V(12)).^2,1));
result_piece4 = max_activity * (1 - (diff_ILD_piece4/V(start_ILDwidth+3)));

%IA_piece5
V_mat = repmat(V(start_ILDpeakshift5:end_ILDpeakshift5)',1,size(ILD_piece2,2));
ILD_opt_piece5 = repmat(ILD_piece5(:,IA_bestloc),1,size(ILD_piece5,2));
diff_ILD_piece5 = sqrt(mean((ILD_piece5 - ILD_opt_piece5 + V_mat).^2,1));
result_piece5 = max_activity * (1 - (diff_ILD_piece5/V(start_ILDwidth+4)));

%IA_piece6
V_mat = repmat(V(start_ILDpeakshift6:end_ILDpeakshift6)',1,size(ILD_piece2,2));
ILD_opt_piece6 = repmat(ILD_piece6(:,IA_bestloc),1,size(ILD_piece6,2));
diff_ILD_piece6 = sqrt(mean((ILD_piece6 - ILD_opt_piece6 + V_mat).^2,1));
result_piece6 = max_activity * (1 - (diff_ILD_piece6/V(start_ILDwidth+5)));

%IA_piece7
V_mat = repmat(V(start_ILDpeakshift7:end_ILDpeakshift7)',1,size(ILD_piece2,2));
ILD_opt_piece7 = repmat(ILD_piece7(:,IA_bestloc),1,size(ILD_piece7,2));
diff_ILD_piece7 = sqrt(mean((ILD_piece7 - ILD_opt_piece7 + V(15)).^2,1));
result_piece7 = max_activity * (1 - (diff_ILD_piece7/V(start_ILDwidth+6)));

%IA_piece8
V_mat = repmat(V(start_ILDpeakshift8:end_ILDpeakshift8)',1,size(ILD_piece2,2));
ILD_opt_piece8 = repmat(ILD_piece8(:,IA_bestloc),1,size(ILD_piece8,2));
diff_ILD_piece8 = sqrt(mean((ILD_piece8 - ILD_opt_piece8 + V_mat).^2,1));
result_piece8 = max_activity * (1 - (diff_ILD_piece8/V(start_ILDwidth+7)));


result = (V(1)*result_piece1) + (V(2)*result_piece2)...
   + (V(3)*result_piece3) + (V(4)*result_piece4) + (V(5)*result_piece5)...
   + (V(6)*result_piece6) + (V(7)*result_piece7) + (V(8)*result_piece8);

index = find(result < 0);
result(index) = 0;

temp = corrcoef(IA_fullspec,result);
cc = temp(1,2);
%err = sum(((result - IA_fullspec).^2))/sum((IA_fullspec).^2);
err = 100 - (100*(cc.^2));

return