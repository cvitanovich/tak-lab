function err = errorfun_spectralintegration2(V,freq_axis,Tonal_meansurf,BP_meansurf)
%Error function for input to fminsearch -- fit the Tonal and BandPass data to the model developed
%by Nelken, Kim and Young (J Neurophysiol 78(1997)800-811)
%var1 is the effective bandwidth for equalization of power
%var2 is the exponential term for the Tonal response, prior to summing

target_bandwidth = 1/3;

%determine weights
f_weights(1) = 1; %from paper
f_weights(length(freq_axis)) = 1;
for k = 2:(length(freq_axis) - 1)
   f_weights(k) = (freq_axis(k+1) - freq_axis(k-1)) / 2;
end
f_weights = repmat(f_weights',1,size(Tonal_meansurf,2));


for num_freq = 1:length(freq_axis)
   [lo_limit,hi_limit] = bandlimits(freq_axis(num_freq),target_bandwidth);
   [temp,index_1] = min(abs(freq_axis - lo_limit));
   [temp,index_2] = min(abs(freq_axis - hi_limit));
   
   %calculate the new ILD/Freq surface from the Tonal data
   Tonal_meansurf_new(num_freq,:) = (1/V(1)) .* ...
      sum(f_weights(index_1:index_2,:) .* ...
      ((Tonal_meansurf(index_1:index_2,:)).^1),1);
end

cc = corrcoef(Tonal_meansurf_new(:),BP_meansurf(:));
err = 1 - cc(1,2);

