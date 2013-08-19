function [ERBforward,ERBfeedback,response,freqScale] = oct0897(fs,numChannels,fftpts)

% [ERBforward,ERBfeedback,response,freqScale] = oct0897(fs,numChannels,fftpts)
% Oct0897
% script to run ERBfiltering and check it out 
% with an impulse response
% defaults:
%      fftpts = 1024*8
%      numChannels = 0 forces ERBfilt3 to use 1/32 octave filtering
%      fs = 30000

if (nargin < 3)   fftpts=1024*8;    end;
if (nargin < 2)   numChannels = 0;  end;           % = 0 forces ERBfilt3 to use 1/32 octave spacing
if (nargin < 1)   fs = 30000;       end;

lowFreq = 1000;     
impulse = [1 zeros(1,1023)];
[ERBforward,ERBfeedback,cf] = ERBfilt3(fs,numChannels,lowFreq);
y = filtbank(ERBforward,ERBfeedback,impulse);
figure
subplot(2,1,1)
response = 20*log10(abs(fft(y',fftpts)));
freqScale = (0:fftpts-1)/fftpts*fs;
axis([2 4 -70 10]);
semilogx(freqScale(1:fftpts/2),response(1:fftpts/2,:));
ylabel ('dB magnitude')
xlabel ('Frequency Hz')
title ('ERB gammatone filter bank')

% compare Q10dB pts
[m,n] = size(response);
for i = 1:n
  index = find(response(1:fftpts/2,i)>=-10.5);
  Q10(i)=cf(i)/(max(freqScale(index))-min(freqScale(index)));
end
subplot(2,1,2)
semilogx(cf,Q10)
hold on
semilogx(cf,.074*(cf.^.504),'r');    % Koppl's
grid on
ylabel ('Q10 dB')
xlabel ('Frequency Hz')
title ('Q10 compared Koppl (red) vs filterbank (blue)')

% compare Q40dB pts

%figure
%for i = 1:20
%  index = find(response(1:fftpts/2,i)>=-40.5);
%  Q40(i)=cf(i)/(max(freqScale(index))-min(freqScale(index)));
%end

%semilogx(cf,Q40)
%hold on
%semilogx(cf,.036*(cf.^.502),'r');   % Koppl's
%grid on
