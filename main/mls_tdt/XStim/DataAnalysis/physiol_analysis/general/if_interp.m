function [new_if_surface,new_if_freqaxis,ILDmat_freqind] = ...
   if_interp(if_surface,if_freqaxis,if_ildaxis,ILDmat_freqaxis)

maxfreq = max(if_freqaxis);
minfreq = min(if_freqaxis);

if(min(ILDmat_freqaxis) > minfreq | max(ILDmat_freqaxis) < maxfreq)
   disp('Frequency band of HRTFs is not wide enought for specified');
   disp(['frequency limits fmin = ' num2str(minfreq) ...
         ' and fmax = ' num2str(maxfreq)])
   return
end

[y,temp1] = min(abs(ILDmat_freqaxis - minfreq));
[y,temp2] = min(abs(ILDmat_freqaxis - maxfreq));
ILDmat_freqind = [temp1 temp2];

new_if_freqaxis = ILDmat_freqaxis(ILDmat_freqind(1):ILDmat_freqind(2));
[x,y] = meshgrid(if_freqaxis,if_ildaxis');
[xi,yi] = meshgrid(new_if_freqaxis,if_ildaxis);
new_if_surface = ...
   interp2(x,y,if_surface',xi,yi,'spline');
new_if_surface = new_if_surface';

return
