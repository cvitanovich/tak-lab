function [Neuron] = iftransform_4(Neuron,ILDspect_opt,HRTFinfo,measurement_type,ILDspect_alt)
%Transforms an ILD-Freq RS taken with a spatial pre-filter into
%a true ILD-Freq RS

%Set stringency for identifying differences from optimal ILD spectrum spike rate
sigma_criterion = 0.5;

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

%Get new ild axis if necessary
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

switch measurement_type
case {'gsif'} %Gammatone Stack ILD-Freq
   Neuron.gsif_ildaxis = newildaxis;
   Neuron.gsif_raw_meansurf = tempmeansurf;
   Neuron.gsif_orig_meansurf = meansurf - zeroild_mean;
   Neuron.gsif_orig_stdsurf = stdsurf;
case {'bbif'} %BroadBand Noise ILD-Freq
   Neuron.bbif_ildaxis = newildaxis;
   Neuron.bbif_raw_meansurf = tempmeansurf;
   Neuron.bbif_orig_meansurf = meansurf - zeroild_mean;
   Neuron.bbif_orig_stdsurf = stdsurf;
end

tempmeansurf = meansurf;
[Xnew,Ynew] = meshgrid(freqaxis,newildaxis);
negsurf = meansurf - zeroild_mean;
neg_index = find(negsurf < (-1.*sigma_criterion.*stdsurf));
out_index = find(negsurf >= (-1.*sigma_criterion.*stdsurf) & (Ynew' < min(ildaxis) | Ynew' > max(ildaxis)));
positive_index = find(negsurf >= (-1.*sigma_criterion.*stdsurf) & (Ynew' >= min(ildaxis) & Ynew' <= max(ildaxis)));
tempmeansurf(neg_index) = negsurf(neg_index);
tempmeansurf(out_index) = negsurf(out_index);
%positive_index = find(negsurf >= (-1.*sigma_criterion.*stdsurf) & ~isnan(tempmeansurf));
%nan_index = find(isnan(tempmeansurf));

%Get Z-scores
Zscores = negsurf./stdsurf;

if(nargin > 4)
   ILDspect_opt = ILDspect_alt;
end




%Fourth, apply the multiplication filter to the appropriate tonal or gammatonal RS
if(isfield(Neuron,'tif_mean') & ...
      ~isfield(Neuron,'gif_mean') & ...
      (strcmp(measurement_type,'gif') | strcmp(measurement_type,'tsif') | ...
      strcmp(measurement_type,'gsif') | strcmp(measurement_type,'bbif')))
   %Interpolate mean to match indices
   [Xold,Yold] = meshgrid(Neuron.tif_freqaxis,Neuron.tif_ildaxis);
   newmeansurf = interp2(Yold',Xold',Neuron.tif_mean,Ynew',Xnew');
   tempmeansurf(positive_index) = tempmeansurf(positive_index) .* newmeansurf(positive_index);
%   tempmeansurf(nan_index) = Neuron.tif_mean(nan_index);
end
if(isfield(Neuron,'gif_mean') & (strcmp(measurement_type,'gsif') | strcmp(measurement_type,'tsif') | strcmp(measurement_type,'bbif')))
   [Xold,Yold] = meshgrid(Neuron.gif_freqaxis,Neuron.gif_ildaxis);
   newmeansurf = interp2(Yold',Xold',Neuron.gif_mean,Ynew',Xnew');
   tempmeansurf(positive_index) = tempmeansurf(positive_index) .* newmeansurf(positive_index);
%   tempmeansurf(nan_index) = Neuron.gif_mean(nan_index);
end

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