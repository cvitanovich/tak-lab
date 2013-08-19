% script for ramped and damped sinusoids 
% from RD Patterson
% JASA 96(3):1412

power2 = 2;        % determines halflife in msec and amplitude
c = -.693147;      % scaling constant
T = 25;            % rep period in msec
w = 2*pi*1000;     % freq in radians
hl = 2^power2;     % half life in msec
A = 1 + power2*sqrt(2);
                   % amplitude needs to be increased
                   % by sqrt(2)
                   % each time the half-life is doubled

% for one damped or ramped rep:
i = 0;
time = 0:.1:T-.1;    % time in msec
for t = time
   i = i+1;
   damp(i) = A * exp(c*t/hl)*sin(w*t/1000);
   ramp(i) = A * exp(c*(T-t)/hl)*sin(w*(T-t)/1000);
end

% to see two cycles:

damp = [damp damp];
ramp = [ramp ramp];
totaltime = 0:.1:(2*T)-.1;
plot(totaltime,damp);
hold on
plot(totaltime,ramp,'r');

% calculate ERB filter bank coefficients
% 10000 Hz sampling rate
% 80 filters
% lowest freq = 100 Hz

[forward, feedback]=ERBfilt(10000,80,100);

% do the filtering
yramp = filtbank(forward,feedback,ramp);
ydamp = filtbank(forward,feedback,damp);

% plot damped version
figure
for i = 1:80
  plot(ydamp(i,:)-(i*.1*power2));
  hold on
end

% plot ramped version
figure
for i = 1:80
  plot(yramp(i,:)-(i*power2*.1));
  hold on
end

% plot one of each at freq slightly above CF
figure
plot(yramp(10,:))
hold on
plot(ydamp(10,:),'r')

% plot peak-to-peak values
figure
plot(abs(yramp(10,:)))
hold on
plot(abs(ydamp(10,:)),'r')


% to test filterbank on noise
 numChannels = 0;
 lowFreq = 1000;
 fs = 30000;
 fftpts=1024;

 [forward, feedback, cf]=ERBfilt3(fs,numChannels,lowFreq);
 x = r_noi('e:\kip\noisetc\t30\forward.t30');
 y = filtbank(forward,feedback,x);
 % channels 4 to 15 are of interest here
 figure
response = 20*log10(abs(fft(y',fftpts)));
freqScale = (0:fftpts-1)/fftpts*fs;
axis([2 4 -70 10]);
semilogx(freqScale(1:fftpts/2),response(1:fftpts/2,:));
 
 

 % to test ERB filter bank on impulse:
 numChannels = 0;
 lowFreq = 1000;
 fs = 30000;
 fftpts=1024;

impulse = [1 zeros(1,1023)];
[ERBforward,ERBfeedback,cf] = ERBfilt3(fs,numChannels,lowFreq);
y = filtbank(ERBforward,ERBfeedback,impulse);
figure
response = 20*log10(abs(fft(y',fftpts)));
freqScale = (0:fftpts-1)/fftpts*fs;
axis([2 4 -70 10]);
semilogx(freqScale(1:fftpts/2),response(1:fftpts/2,:));
