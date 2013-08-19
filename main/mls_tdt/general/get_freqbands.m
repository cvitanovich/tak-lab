function [freq_bands] = get_freqbands(limits,bandwidth)
% Function to provide incremental frequency bands of a specified bandwidth for specified broadband limits
% limits: [lo hi], in Hz
% bandwidth: in octaves; e.g., 1/3, 1/8, 1/2

freq_bands = round([limits(1) limits(1)*(2^bandwidth)]);
band_num = 1;
while freq_bands(band_num,2) < limits(2)
   band_num = band_num + 1;
   freq_bands(band_num,1) = round(freq_bands(band_num-1,2));
   freq_bands(band_num,2) = round(freq_bands(band_num,1) * (2^bandwidth));
end
if(freq_bands(band_num,2) > limits(2)) freq_bands(band_num,2) = round(limits(2)); end
if(freq_bands(band_num,2) < freq_bands(band_num,1))
   disp('Error: limits must be specified in increasing order.')
   freq_bands = [];
end
