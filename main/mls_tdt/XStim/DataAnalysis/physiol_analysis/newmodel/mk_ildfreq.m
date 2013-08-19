function [result,V,freq_weight] = ...
   mk_ildfreq(V,...
   if_surface,...
   if_freqaxis,...
   if_ildaxis,...
   ILD_matrix,...
   HRTF_freqaxis,...
   IA_bestloc_ind,...
   freq_bands)

if(isempty(freq_bands))
   freq_bands(:,1:2) = [if_freqaxis if_freqaxis];
end

%specify indexing for V
start_ILDwidth = 1; end_ILDwidth = start_ILDwidth + size(freq_bands,1) - 1;
start_ILDshift = end_ILDwidth + 1; end_ILDshift = start_ILDshift + size(freq_bands,1) - 1;
lenV = length(V);

%Assume that ILDwidths must be less than a certain dB and > 0
max_allowable_ILDwidth = 25; %dB
temp = V(start_ILDwidth:end_ILDwidth);
index_ILDwidth = find(temp > max_allowable_ILDwidth);
if(~isempty(index_ILDwidth))
   temp(index_ILDwidth) = max_allowable_ILDwidth;
   V(start_ILDwidth:end_ILDwidth) = temp;
end
temp = V(start_ILDwidth:end_ILDwidth);
index_ILDwidth = find(temp <= 0);
if(~isempty(index_ILDwidth))
   temp(index_ILDwidth) = 1;
   V(start_ILDwidth:end_ILDwidth) = temp;
end

%Assume that ILDopt_shift must be less than  +/- some amount
max_allowable_ILDshift = 5; %dB
temp = V(start_ILDshift:end_ILDshift);
index_ILDshift = find(abs(temp) > max_allowable_ILDshift);
if(~isempty(index_ILDshift))
   temp(index_ILDshift) = sign(temp(index_ILDshift)) * max_allowable_ILDshift;
   V(start_ILDshift:end_ILDshift) = temp;
end

for band_num = 1:size(freq_bands,1)
   [y,minfreqindex1] = min(abs(HRTF_freqaxis - freq_bands(band_num,1)));
   [y,maxfreqindex1] = min(abs(HRTF_freqaxis - freq_bands(band_num,2)));
   [y,minfreqindex2] = min(abs(if_freqaxis - freq_bands(band_num,1)));
   [y,maxfreqindex2] = min(abs(if_freqaxis - freq_bands(band_num,2)));
   
   freq_weight(band_num) = mean(max(if_surface(minfreqindex2:maxfreqindex2,:),[],1));
   
   ILD_opt = repmat(ILD_matrix(minfreqindex1:maxfreqindex1,IA_bestloc_ind),1,size(ILD_matrix,2));
   dist = sqrt(mean((ILD_matrix(minfreqindex1:maxfreqindex1,:) -...
      ILD_opt + V(start_ILDshift - 1 + band_num)).^2,1));
   [ILD_matrix_freqband(band_num,:),index_ILDsort] = ...
      sort(mean(ILD_matrix(minfreqindex1:maxfreqindex1,:),1));
   clear temp temp2
   temp = freq_weight(band_num) * ...
      (1 - ((dist/V(start_ILDwidth-1+band_num)).^V(lenV)));
   temp2 = interp1(ILD_matrix_freqband(band_num,:),...
      temp(index_ILDsort),if_ildaxis);
   result_pieces(band_num,:) = temp2';
   index_nan = find(isnan(result_pieces(band_num,:)));
   result_pieces(band_num,index_nan) = 0;
end

freq_step = freq_bands(2,1) - freq_bands(1,1);
result = interp2((freq_bands(1,1):...
   freq_step:...
   freq_bands(size(freq_bands,1),1)),...
   if_ildaxis,result_pieces',if_freqaxis,if_ildaxis');

index = find(result < 0);
result(index) = 0;

clear index
index = find(isnan(result));
result(index) = 0;


return
