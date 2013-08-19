function [diam_itd, az, el, wa_itd] = itd2space3(itd_data, itd_ax,...
hrtfitd, hrtflocs, hrtffreqs, freq_curve, freq_axis);
% [diam_itd, az, el, wa_itd] = itd2space3(itd_data, itd_ax, hrtfitd,
%             hrtflocs, hrtffreqs, freq_curve, freq_axis);
% ITD2SPACE3  convert itd tuning curve into spatial receptive field
%             using hrtf-based ITD values.  This version does this by
%             using the frequency-specific itd values and doing a
%             frequency-weighted
%             average to find the best estimate of 'effective itd' at each
%             location.
%             hrtfitd is an array of frequency-specific ITD values,
%             locations along rows, freqs along columns
%             hrtflocs and hrtffreqs are the axes for hrtfitd
%             itd_data is the ITD tuning curve, and itd_ax is the
%             corresponding axes values
%             freq_curve is the frequency response curve (single row
%             array) and axis for that curve.
%             output wa_itd is the weighted average ITD surface

% make the freq curve and itd spectrum cover the same frequencies
itdminfreq = max(min(hrtffreqs),min(freq_axis));
itdmaxfreq = min(max(hrtffreqs),max(freq_axis));
redi = hrtffreqs<=itdmaxfreq & hrtffreqs>=itdminfreq;
red_hrtfitd = hrtfitd(:,redi);
red_hrtfitd_freq_axis = hrtffreqs(redi);
new_freq_curve = interp1(freq_axis, freq_curve, red_hrtfitd_freq_axis);

% weighted sum
freq_curve2d = (ones(size(red_hrtfitd,1),1))*new_freq_curve;
witd = freq_curve2d.*red_hrtfitd;
itdavg = sum(witd')/sum(new_freq_curve); 
itdavg = itdavg';

% interpolate
space_itd = interp1(itd_ax, itd_data, itdavg);

%plot(itd_ax, itd_data);
%hold on;
%plot(hrtfitd,space_itd,'*r');
%hold off;

% set all out of range values to the minimum value in the itd curve
% not a perfect solution, i know
space_itd(isnan(space_itd)) = min(space_itd)*ones(size(space_itd(isnan(space_itd))));

% form into diamond
[az,el,diam_itd] = array2diamond(space_itd, hrtflocs);    

if (nargout>3)
   [az,el,wa_itd] = array2diamond(itdavg, hrtflocs);

end;

return; 