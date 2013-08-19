cf=1825;
fs= 30000;
freq = 0:15000/1023:15000;
impulse = [1 zeros(1,1023)];
mini=8;
maxi = 683;


oct1097

figure
[forward, feedback, cf]=ERBfilt2(fs,cf);
 y = filter(forward,feedback,impulse);
 yabs = abs(rfft(y,1024*2));
 mindB = min(HIdb);
 semilogx(freq(mini:maxi), -(20*log10(yabs(mini:maxi))- mindB))
 hold on
 semilogx(HIfreq, HIdb,'r.')
 
 cf=615;
 [forward, feedback, cf]=ERBfilt2(fs,cf);
 y = filter(forward,feedback,impulse);
 yabs = abs(rfft(y,1024*2));
 mindB = min(LOdb);
 semilogx(freq(mini:maxi), -(20*log10(yabs(mini:maxi))- mindB))
 semilogx(LOfreq, LOdb,'r.')
 
 
 cf=5000;
 [forward, feedback, cf]=ERBfilt2(fs,cf);
 y = filter(forward,feedback,impulse);
 yabs = abs(rfft(y,1024*2));
 mindB = min(TOPdb);
 semilogx(freq(mini:maxi), -(20*log10(yabs(mini:maxi))- mindB))
 semilogx(TOPfreq, TOPdb,'r.')
