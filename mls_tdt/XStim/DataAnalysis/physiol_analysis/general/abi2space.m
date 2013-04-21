function [ABI_diamond, abi_az, abi_el] = abi2space(ABI_matrix, ABI_freqs, ABI_locations, abi_axis,...
   abi_curve, Tonal_freqcurve, Tonal_freqaxis, cell_threshhold)

%Convert the ABI matrix to values that take into account the actual dB level of the sound
ABI_matrix = ABI_matrix + (cell_threshhold + 20); %sounds were played at 20 dB over threshhold

%Make it so the Tonal frequency curve and the ABI matrix cover the same frequencies
minfreq = max(min(ABI_freqs),min(Tonal_freqaxis));
maxfreq = min(max(ABI_freqs),max(Tonal_freqaxis));
reduced_ABIfindex = ABI_freqs >= minfreq & ABI_freqs <= maxfreq;
reduced_ABI_matrix = ABI_matrix(reduced_ABIfindex,:);
reduced_ABI_freqs = ABI_freqs(reduced_ABIfindex);
new_Tonal_freqcurve = interp1(Tonal_freqaxis,Tonal_freqcurve,reduced_ABI_freqs);

%Perform a weighted summation in order to obtain the average ABI responses
Tonal_freqcurve_matrix = new_Tonal_freqcurve*(ones(1,size(reduced_ABI_matrix,2)));
Weighted_ABI_matrix = Tonal_freqcurve_matrix .* reduced_ABI_matrix;
ABI_avgmatrix = sum(Weighted_ABI_matrix,1)/sum(new_Tonal_freqcurve);

%Perform interpolation to get the abi response vector
spatial_ABIresponse = interp1(abi_axis, abi_curve, ABI_avgmatrix);

%Set all out-of-range values to minimum value in abi_curve
spatial_ABIresponse(isnan(spatial_ABIresponse)) =...
   min(spatial_ABIresponse)*ones(size(spatial_ABIresponse(isnan(spatial_ABIresponse))));

%Convert to diamond
[abi_az,abi_el,ABI_diamond] = array2diamond(spatial_ABIresponse,ABI_locations);
ABI_diamond = ABI_diamond./max(max(abs(ABI_diamond)));

return