function [Neuron] = iftransform(Neuron,ILDspect_opt,HRTFinfo,measurement_type,ILDspect_alt)
%Transforms an ILD-Freq RS taken with a spatial pre-filter into
%a true ILD-Freq RS

%Set stringency for identifying differences from optimal ILD spectrum spike rate
sigma_criterion = 1;
num_sigdiff_criterion = 4;

switch measurement_type
case {'tif'} %Tonal ILD-Freq
   freqaxis  = Neuron.tif_freqaxis;
   ildaxis   = Neuron.tif_ildaxis;
   meansurf  = Neuron.tif_mean;
   stdsurf   = Neuron.tif_std;
case {'tsif'} %Tone Stack ILD-Freq
   freqaxis  = Neuron.tsif_freqaxis;
   ildaxis   = Neuron.tsif_ildaxis;
   meansurf  = Neuron.tsif_mean;
   stdsurf   = Neuron.tsif_std;
case {'gif'} %Gammatone ILD-Freq
   freqaxis  = Neuron.gif_freqaxis;
   ildaxis   = Neuron.gif_ildaxis;
   meansurf  = Neuron.gif_mean;
   stdsurf   = Neuron.gif_std;
case {'gsif'} %Gammatone Stack ILD-Freq
   freqaxis  = Neuron.gsif_freqaxis;
   ildaxis   = Neuron.gsif_ildaxis;
   meansurf  = Neuron.gsif_mean;
   stdsurf   = Neuron.gsif_std;
otherwise
   error('Unknown measurement type.')
end

%First, transform by correcting to true ILD
tempmeansurf = meansurf;
%Get the ild index for 0 ild
zeroildind = find(ildaxis == 0);
tempmeansurf = tempmeansurf - mean(tempmeansurf(:,zeroildind));
tempstdsurf  = stdsurf;

clear meansurf stdsurf
meansurf = NaN * ones(size(tempmeansurf));
stdsurf  = NaN * ones(size(tempstdsurf));

for freqnum = 1:length(freqaxis)
   %Get the frequency index for the ILD spectrum
   [y,HRTFfreqind] = min(abs(HRTFinfo.hrtf_freqs - freqaxis(freqnum)));
   for ildnum = 1:length(ildaxis)
      trueild(freqnum,ildnum) = ildaxis(ildnum) - ILDspect_opt(HRTFfreqind);
      [y,ildind] = min(abs(ildaxis - trueild(freqnum,ildnum)))
      if(y > 5)
         meansurf(freqnum,ildind) = NaN;
      else
         meansurf(freqnum,ildind) = tempmeansurf(freqnum,ildnum);
         stdsurf(freqnum,ildind)  = tempstdsurf(freqnum,ildnum);
      end
   end
end


switch measurement_type
case {'tif'} %Tonal ILD-Freq
   Neuron.tif_mean = meansurf;
   Neuron.tif_std = stdsurf;
case {'tsif'} %Tone Stack ILD-Freq
   Neuron.tsif_mean = meansurf;
   Neuron.tsif_std = stdsurf;
case {'gif'} %Gammatone ILD-Freq
   Neuron.gif_mean = meansurf;
   Neuron.gsif_std = stdsurf;
case {'gsif'} %Gammatone Stack ILD-Freq
   Neuron.gsif_mean = meansurf;
   Neuron.gsif_std = stdsurf;
otherwise
   error('Unknown measurement type.')
end

return