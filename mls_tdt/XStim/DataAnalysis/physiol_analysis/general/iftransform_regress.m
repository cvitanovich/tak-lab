function [new_bbif_meansurf,newildaxis] = ...
   iftransform_regress(bbif_meansurf,freqaxis,ildaxis,ILDspect_opt,HRTFinfo)
%Transforms an ILD-Freq RS taken with a spatial pre-filter into
%a true ILD-Freq RS



%First, transform by correcting to true ILD
meansurf = bbif_meansurf;
tempmeansurf = meansurf;
clear meansurf

%Expand ild axis if necessary
maxild = max(ildaxis); minild = min(ildaxis);
for freq_num = 1:length(freqaxis)
   %Get the frequency indices for the ILD spectrum
   [freq_bands(freq_num,1),freq_bands(freq_num,2)] = bandlimits(freqaxis(freq_num),1/12);
         if(freq_bands(freq_num,1) < 2000) freq_bands(freq_num,1) = 2000; end
         if(freq_bands(freq_num,2) < 2000) freq_bands(freq_num,2) = 2000; end
         if(freq_bands(freq_num,1) > 11000) freq_bands(freq_num,1) = 11000; end
         if(freq_bands(freq_num,2) > 11000) freq_bands(freq_num,2) = 11000; end
   [y,lofreqind] = min(abs(HRTFinfo.hrtf_freqs - freq_bands(freq_num,1)));
   [y,hifreqind] = min(abs(HRTFinfo.hrtf_freqs - freq_bands(freq_num,2)));
   sign_RMS_ILD = sign(mean(ILDspect_opt(lofreqind:hifreqind))) * ...
      rms(ILDspect_opt(lofreqind:hifreqind));
   for ildnum = 1:length(ildaxis)
      trueild = ildaxis(ildnum) + sign_RMS_ILD;
      [y,ildind] = min(abs(ildaxis - trueild));
      if(y > 1 & trueild > maxild)
         maxild = trueild;
      elseif (y > 1 & trueild < minild)
         minild = trueild;
      end
   end
end

newildaxis = round(minild):1:round(maxild);

%Rearrange the ILD-Frequency plot
meansurf = NaN * ones(size(tempmeansurf,1),length(newildaxis));
stdsurf  = NaN * ones(size(tempmeansurf,1),length(newildaxis));

for freq_num = 1:length(freqaxis)
   %Get the frequency indices for the ILD spectrum
   [freq_bands(freq_num,1),freq_bands(freq_num,2)] = bandlimits(freqaxis(freq_num),1/12);
         if(freq_bands(freq_num,1) < 2000) freq_bands(freq_num,1) = 2000; end
         if(freq_bands(freq_num,2) < 2000) freq_bands(freq_num,2) = 2000; end
         if(freq_bands(freq_num,1) > 11000) freq_bands(freq_num,1) = 11000; end
         if(freq_bands(freq_num,2) > 11000) freq_bands(freq_num,2) = 11000; end
   [y,lofreqind] = min(abs(HRTFinfo.hrtf_freqs - freq_bands(freq_num,1)));
   [y,hifreqind] = min(abs(HRTFinfo.hrtf_freqs - freq_bands(freq_num,2)));
   sign_RMS_ILD = sign(mean(ILDspect_opt(lofreqind:hifreqind))) * ...
      rms(ILDspect_opt(lofreqind:hifreqind));
   for ildnum = 1:length(ildaxis)
      trueild = ildaxis(ildnum) + sign_RMS_ILD;
      [y,ildind] = min(abs(newildaxis - trueild));
      if(y > 1)
         meansurf(freq_num,ildind) = NaN;
      else
         meansurf(freq_num,ildind) = tempmeansurf(freq_num,ildnum);
      end
   end
end

new_bbif_meansurf = meansurf;
return