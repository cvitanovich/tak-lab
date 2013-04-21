function [Neuron] = iftransform(Neuron,ILDspect_opt,HRTFinfo,measurement_type,ILDspect_alt)
%Transforms an ILD-Freq RS taken with a spatial pre-filter into
%a true ILD-Freq RS

%Set stringency for identifying differences from optimal ILD spectrum spike rate
sigma_criterion = 1;
num_sigdiff_criterion = 2;

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
tempstdsurf  = stdsurf;
clear meansurf stdsurf
meansurf = NaN * ones(size(tempmeansurf));
stdsurf  = NaN * ones(size(tempstdsurf));

for freqnum = 1:length(freqaxis)
   %Get the frequency index for the ILD spectrum
   [y,HRTFfreqind] = min(abs(HRTFinfo.hrtf_freqs - freqaxis(freqnum)));
   for ildnum = 1:length(ildaxis)
      trueild = ildaxis(ildnum) - ILDspect_opt(HRTFfreqind);
      [y,ildind] = min(abs(ildaxis - trueild));
      if(y > 5)
         meansurf(freqnum,ildind) = NaN;
      else
         meansurf(freqnum,ildind) = tempmeansurf(freqnum,ildnum);
         stdsurf(freqnum,ildind)  = tempstdsurf(freqnum,ildnum);
      end
   end
end


if(nargin > 4)
   ILDspect_opt = ILDspect_alt;
end

%Second, transform by identifying maxima and relating all other ILD results to those
for freqnum = 1:length(freqaxis)
   %Get the frequency index for the ILD spectrum
   [y,HRTFfreqind] = min(abs(HRTFinfo.hrtf_freqs - freqaxis(freqnum)));
   %Get the ildaxis index closest to the optimal ILD value for this frequency band
   [y,ildind_opt]  = min(abs(ildaxis - ILDspect_opt(HRTFfreqind)));
   %Get the maximal activity for current frequency
   [ratemax,ildindmax] = max(meansurf(freqnum,:));
   %Make sure that if there are >1 maxima, take the one closest to the optimal ILD value
   if(length(ratemax) > 1)
      [y,ildind2] = min(abs(ildaxis(ildindmax) - ILDspect_opt(HRTFfreqind)));
      temp1 = ratemax(ildind2);
      clear ratemax; 
      ratemax = temp1;
      temp1 = ildindmax(ildind2);
      clear ildindmax
      ildindmax = temp1;
   end
   maxflag = 1;
   for ildnum = 1:length(ildaxis)
      if(ildnum == ildindmax)
         tempmeansurf(freqnum,ildnum) = 1;
      elseif( abs(meansurf(freqnum,ildindmax) - meansurf(freqnum,ildnum)) >= ...
            (sigma_criterion * stdsurf(freqnum,ildnum)) )
         tempmeansurf(freqnum,ildnum) = ...
            (meansurf(freqnum,ildnum) - meansurf(freqnum,ildindmax))/meansurf(freqnum,ildindmax);
      else
         if(ildnum == ildind_opt)
            tempmeansurf(freqnum,ildnum) = 1;
            maxflag = 0;
         else
            tempmeansurf(freqnum,ildnum) = 0;
         end
      end
   end
   tempmeansurf(freqnum,ildindmax) = maxflag * tempmeansurf(freqnum,ildindmax);
end

%Third, make sure to only weight those frequency bands with measurable changes,
%determined by the statistical criterion set above
max_numsigdiff = -1;
for freqnum = 1:length(freqaxis)
   numsigdiff(freqnum) = length(find(tempmeansurf(freqnum,:) < 0));
   if(numsigdiff(freqnum) < num_sigdiff_criterion)
      tempmeansurf(freqnum,:) = zeros(size(tempmeansurf(freqnum,:)));
   end
end
multmat = repmat(numsigdiff'/max(numsigdiff),1,length(ildaxis));
tempmeansurf = tempmeansurf .* multmat;

switch measurement_type
case {'tif'} %Tonal ILD-Freq
   Neuron.tif_mean = tempmeansurf;
   Neuron.tif_std = stdsurf;
case {'tsif'} %Tone Stack ILD-Freq
   Neuron.tsif_mean = tempmeansurf;
   Neuron.tsif_std = stdsurf;
case {'gif'} %Gammatone ILD-Freq
   Neuron.gif_mean = tempmeansurf;
   Neuron.gsif_std = stdsurf;
case {'gsif'} %Gammatone Stack ILD-Freq
   Neuron.gsif_mean = tempmeansurf;
   Neuron.gsif_std = stdsurf;
otherwise
   error('Unknown measurement type.')
end

return