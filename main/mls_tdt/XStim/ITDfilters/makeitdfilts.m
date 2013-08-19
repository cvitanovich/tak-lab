%MaKeITDfilts
%Make the ITD filters necessary for use with XStim
%ITD range between -250 and 250 us

clear

Fs = 30000;
Fs2 = 1e6;
filt_length = 256;
temp = zeros(1,filt_length);
temp(128) = 1;

itd_min = 0;
itd_max = 250;


for itd = itd_min:itd_max
   clear itd_filt lenupsample itd_filt_up itd_filt_dn itd_filt_temp
   itd_filt = temp;
   itd_filt_up = resample(itd_filt,Fs2,Fs); %upsample
   lenupsample = length(itd_filt_up);
   itd_filt_up = [zeros(1,itd) itd_filt_up];
   itd_filt_temp = resample(itd_filt_up(1:lenupsample),Fs,Fs2); 
   itd_filt_dn = itd_filt_temp(1:filt_length);
   
   itd_filt_dn = itd_filt_dn/max(itd_filt_dn);
   eval(['itd_filt' num2str(itd) ' = itd_filt_dn;']);
   eval(['save ' FN.ITD_path 'itdfilt' num2str(itd) ' itd_filt' num2str(itd)]);
   
   
   disp(['Finished ITD = ' num2str(itd)])
end

figure
plot(0:Fs2/Fs2:lenupsample-(Fs2/Fs2),itd_filt_up(1:lenupsample));
hold on
plot(0:Fs2/Fs:lenupsample - (Fs2/Fs),itd_filt1,'r');