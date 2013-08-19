function [Neuron] = iftransform_5(Neuron,ILDspect_opt,HRTFinfo,measurement_type,ILDspect_alt)
%Transforms an ILD-Freq RS taken with a spatial pre-filter into
%a true ILD-Freq RS

%Set stringency for identifying differences from optimal ILD spectrum spike rate
sigma_criterion = 1.0;

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
%Discard outlying max and min and take mean of 0 ILD shifted results
dismaxind = find(tempmeansurf(:,zeroild_index) == max(tempmeansurf(:,zeroild_index)));
disminind = find(tempmeansurf(:,zeroild_index) == min(tempmeansurf(:,zeroild_index)));
if(length(dismaxind) == 1 & length(disminind) == 1)
   new_zeroild_array = ...
      tempmeansurf(find(1:size(tempmeansurf,1) ~= dismaxind & ...
      1:size(tempmeansurf,1) ~= disminind),zeroild_index);
elseif(length(dismaxind) == 1 & length(disminind) > 1)
   new_zeroild_array = ...
      tempmeansurf(find(1:size(tempmeansurf,1) ~= dismaxind),zeroild_index);
elseif(length(dismaxind) > 1 & length(disminind) == 1)
   new_zeroild_array = ...
      tempmeansurf(find(1:size(tempmeansurf,1) ~= disminind),zeroild_index);
else
   new_zeroild_array = tempmeansurf(:,zeroild_index);
end
zeroild_mean = nanmean(new_zeroild_array);
tempmeansurf(:,zeroild_index) = zeroild_mean;

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

%Rearrange the ILD-Frequency plot
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
negsurf(6,2) = 0;


switch measurement_type
case {'tif'} %Tonal ILD-Freq
   Neuron.tif_ildaxis = newildaxis;
   Neuron.tif_mean = meansurf;
   Neuron.tif_std = stdsurf;
   return;
case {'gif'} %Gammatone ILD-Freq
   Neuron.gif_ildaxis = newildaxis;
   Neuron.gif_mean = meansurf;
   Neuron.gsif_std = stdsurf;
   return;
case {'tsif'} %Gammatone Stack ILD-Freq
   Neuron.tsif_ildaxis = newildaxis;
   Neuron.tsif_raw_ildaxis = ildaxis;
   Neuron.tsif_raw_meansurf = tempmeansurf;
   Neuron.tsif_orig_meansurf = meansurf - zeroild_mean;
   Neuron.tsif_orig_stdsurf = stdsurf;
case {'gsif'} %Gammatone Stack ILD-Freq
   Neuron.gsif_ildaxis = newildaxis;
   Neuron.gsif_raw_ildaxis = ildaxis;
   Neuron.gsif_raw_meansurf = tempmeansurf;
   Neuron.gsif_orig_meansurf = meansurf - zeroild_mean;
   Neuron.gsif_orig_stdsurf = stdsurf;
case {'bbif'} %BroadBand Noise ILD-Freq
   Neuron.bbif_ildaxis = newildaxis;
   Neuron.bbif_raw_ildaxis = ildaxis;
   Neuron.bbif_raw_meansurf = tempmeansurf;
   Neuron.bbif_orig_meansurf = meansurf - zeroild_mean;
   Neuron.bbif_orig_stdsurf = stdsurf;
end


%Get Z-scores
Zscores = negsurf./stdsurf;

if(nargin > 4)
   ILDspect_opt = ILDspect_alt;
end

%Only do the following for Stacks and Noise
switch measurement_type
case {'tsif','gsif','bbif'}
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
   %Z_diff = abs(negsurf(freqnum,ildindmax) - negsurf(freqnum,ildind_opt))/stdsurf(freqnum,ildind_opt);
   %if(Z_diff < sigma_criterion) ildindmax = ildind_opt; end
   ildindsave(freqnum) = ildindmax;
   for ildnum = 1:length(newildaxis)
      if(ildnum == ildindmax)
         tempmeansurf(freqnum,ildnum) = ratemax;
      else
         if(negsurf(freqnum,ildnum) < 0 & abs(Zscores(freqnum,ildnum)) >= sigma_criterion)
            tempmeansurf(freqnum,ildnum) = negsurf(freqnum,ildnum);
         else
            tempmeansurf(freqnum,ildnum) = 0;
            %tempmeansurf(freqnum,ildnum) = negsurf(freqnum,ildnum);
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
end %end switch

switch measurement_type
case {'tsif'} %Tone Stack ILD-Freq
   Neuron.tsif_mean = tempmeansurf;
   Neuron.tsif_std = stdsurf;
   Neuron.tsif_Zscores = Zscores;
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