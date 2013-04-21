function err = ...
   errorfun_spectralintegration8(V,IA_fullspec,IA_piece1,IA_piece2,ILD_piece1,ILD_piece2,...
   IA_bestloc1,IA_bestloc2)
%This function takes only ILDAlone mean arrays
%errorfun_spectralintegration7(V,IA_fullspec,IA_piece1,IA_piece2,ILD_piece1,ILD_piece2,IA_bestloc)
%V:               variables to fit
%IA_fullspec:     measured ILDAlone RS with broadband noise
%IA_pieces:       bandwidth-limited measured ILDAlone RS's
%ILD_pieces:      ILD matrix from HRTF data (freqs assigned to rows, locs assigned to cols)

ILD_range = min([min(ILD_piece1) min(ILD_piece2)]):max([max(ILD_piece1) max(ILD_piece2)]);

%IA_piece1
for num_freq = 1:size(ILD_piece1,1)
   gauss_vec(num_freq,:) = gaussd(ILD_range,IA_bestloc1(num_freq),V(1));
   for num_loc = 1:size(ILD_piece1,2)
      [y,index] = min(abs(ILD_range - ILD_piece1(num_freq,num_loc)));
      temp(num_freq,num_loc) = gauss_vec(num_freq,index);
   end
end
weights = mean(temp,1);
result_piece1 = IA_piece1 .* weights;

%IA_piece2
for num_freq = 1:size(ILD_piece2,1)
   gauss_vec(num_freq,:) = gaussd(ILD_range,IA_bestloc2(num_freq),V(2));
   for num_loc = 1:size(ILD_piece2,2)
      [y,index] = min(abs(ILD_range - ILD_piece2(num_freq,num_loc)));
      temp(num_freq,num_loc) = gauss_vec(num_freq,index);
   end
end
weights = mean(temp,1);
result_piece2 = IA_piece2 .* weights;

result = result_piece1 + result_piece2;
err = mean( (result - IA_fullspec).^2);

return