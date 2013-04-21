function [Neuron] = iftransform_5(Neuron,ILDspect_opt,HRTFinfo,measurement_type,ILDspect_alt)
%Transforms an ILD-Freq RS taken with a spatial pre-filter into
%a true ILD-Freq RS

%Set stringency for identifying differences from optimal ILD spectrum spike rate
sigma_criterion = 1.6;

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
case {'bbif'} %BroadBand Noise ILD-Freq
   freqaxis  = Neuron.bbif_freqaxis;
   ildaxis   = Neuron.bbif_ildaxis;
   meansurf  = Neuron.bbif_mean;
   stdsurf   = Neuron.bbif_std;
otherwise
   error('Unknown measurement type.')
end

%First, transform by correcting to true ILD
tempmeansurf = meansurf;
tempstdsurf  = stdsurf;
clear meansurf stdsurf

zeroild_index = find(ildaxis == 0);
zeroild_mean = nanmean(tempmeansurf(:,zeroild_index));

%Expand ild axis if necessary
maxild = max(ildaxis); minild = min(ildaxis);
for freqnum = 1:length(freqaxis)
   %Get the frequency index for the ILD spectrum
   [y,HRTFfreqind] = min(abs(HRTFinfo.hrtf_freqs - freqaxis(freqnum)));
   for ildnum = 1:length(ildaxis)
      trueild = ildaxis(ildnum) + ILDspect_opt(HRTFfreqind);
      [y,ildind] = min(abs(ildaxis - trueild));
      if(y > 5 & trueild > maxild)
         maxild = trueild;
      elseif (y > 5 & trueild < minild)
         minild = trueild;
      end
   end
end

newildaxis = round(minild):5:round(maxild);
meansurf = NaN * ones(size(tempmeansurf,1),length(newildaxis));
stdsurf  = NaN * ones(size(tempmeansurf,1),length(newildaxis));


for freqnum = 1:length(freqaxis)
   %Get the frequency index for the ILD spectrum
   [y,HRTFfreqind] = min(abs(HRTFinfo.hrtf_freqs - freqaxis(freqnum)));
   for ildnum = 1:length(ildaxis)
      trueild = ildaxis(ildnum) + ILDspect_opt(HRTFfreqind);
      [y,ildind] = min(abs(newildaxis - trueild));
      if(y > 5)
         meansurf(freqnum,ildind) = NaN;
      else
         meansurf(freqnum,ildind) = tempmeansurf(freqnum,ildnum);
         stdsurf(freqnum,ildind)  = tempstdsurf(freqnum,ildnum);
      end
   end
end

negsurf = meansurf - zeroild_mean;
negsurf(7,11) = 0;


switch measurement_type
case {'gsif'} %Gammatone Stack ILD-Freq
   Neuron.gsif_ildaxis = newildaxis;
   Neuron.gsif_orig_meansurf = negsurf;
   Neuron.gsif_orig_stdsurf = stdsurf;
case {'bbif'} %BroadBand Noise ILD-Freq
   Neuron.bbif_ildaxis = newildaxis;
   Neuron.bbif_orig_meansurf = negsurf;
   Neuron.bbif_orig_stdsurf = stdsurf;
end


%Get Z-scores
Zscores = negsurf./stdsurf;

if(nargin > 4)
   ILDspect_opt = ILDspect_alt;
end

clear tempmeansurf
%Second, transform by identifying maxima and relating all other ILD results to those
for freqnum = 1:length(freqaxis)
   %Get the frequency index for the ILD spectrum
   [y,HRTFfreqind] = min(abs(HRTFinfo.hrtf_freqs - freqaxis(freqnum)));
   %Get the newildaxis index closest to the optimal ILD value for this frequency band
   [y,ildind_opt]  = min(abs(newildaxis - ILDspect_opt(HRTFfreqind)));
   %Get the maximal rate for current frequency
   [ratemax,ildindmax] = max(negsurf(freqnum,:));
   %Make sure that if there are >1 maxima, take the one closest to the optimal ILD value
   if(length(ratemax) > 1)
      [y,ildind2] = min(abs(newildaxis(ildindmax) - ILDspect_opt(HRTFfreqind)));
      temp1 = ratemax(ildind2);
      clear ratemax; 
      ratemax = temp1;
      temp1 = ildindmax(ildind2);
      clear ildindmax
      ildindmax = temp1;
   end
   %Test to see if the maximal rate is significantly different from the rate at the optimal ILD value
   Z_diff = abs(negsurf(freqnum,ildindmax) - negsurf(freqnum,ildind_opt))/stdsurf(freqnum,ildind_opt);
   if(Z_diff < sigma_criterion) ildindmax = ildind_opt; end
   ildindsave(freqnum) = ildindmax;
   for ildnum = 1:length(newildaxis)
      if(ildnum == ildindmax)
         tempmeansurf(freqnum,ildnum) = ratemax;
      else
         if(negsurf(freqnum,ildnum) < 0)
            tempmeansurf(freqnum,ildnum) = negsurf(freqnum,ildnum);
         else
            %tempmeansurf(freqnum,ildnum) = 0;
            tempmeansurf(freqnum,ildnum) = negsurf(freqnum,ildnum);
         end
      end
   end
end

%Third, make sure to weight the frequency bands by the fraction of maximal difference
%across all frequency bands
for freqnum = 1:length(freqaxis)
      if(any(abs(Zscores(freqnum,:)) > sigma_criterion))
         multmat(freqnum,:) = ones(1,length(newildaxis));
      else
         multmat(freqnum,:) = zeros(1,length(newildaxis));
      end
   end
end
tempmeansurf = tempmeansurf .* multmat;

switch measurement_type
case {'tif'} %Tonal ILD-Freq
   Neuron.tif_mean = tempmeansurf;
   Neuron.tif_std = stdsurf;
case {'tsif'} %Tone Stack ILD-Freq
   Neuron.tsif_mean = tempmeansurf;
   Neuron.tsif_std = stdsurf;
   Neuron.tsif_Zscores = Zscores;
case {'gif'} %Gammatone ILD-Freq
   Neuron.gif_mean = tempmeansurf;
   Neuron.gsif_std = stdsurf;
case {'gsif'} %Gammatone Stack ILD-Freq
   Neuron.gsif_mean = tempmeansurf;
   Neuron.gsif_std = stdsurf;
   Neuron.gsif_Zscores = Zscores;
case {'bbif'} %BroadBand Noise ILD-Freq
   Neuron.bbif_mean = tempmeansurf;
   Neuron.bbif_std = stdsurf;
   Neuron.bbif_Zscores = Zscores;
otherwise
   error('Unknown measurement type.')
end

return