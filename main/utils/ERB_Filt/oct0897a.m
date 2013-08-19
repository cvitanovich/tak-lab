% OCT0897a
% script to filter a random noise
% by the Gammatone filterbank
% check out the cross correlations
% and save sounds to files

Fs = 30000;
duration = 0.1;
rampdur = .005;
maxval = 2048;
fpath = 'd:\kip\noisetc\gamma\';
fsuffix = 'b.noi';

% use random noise or some file from disk
Path = input('Enter path for noise or <CR> for internally generated RAND  ','s');
FName = input('Enter filename or (<CR> for internally generated RAND)  ','s');
if isempty(FName)
   noise = rand(1,length(time));
else
   noise = r_noi(Path,FName);
   noise = noise - mean(noise);
   noise = noise/ max([max(noise) abs(min(noise))]);
   duration = length(noise)/Fs;
end
time = 1/Fs:1/Fs:duration;


numChannels = 0;             %this forces ERBfilt3 to use 1/30 octave spacing
% originally 1/6 octave spacing was used for e.g. FIXD
% and 1/32 octave spacing was used for .BPF

lowFreq = 1000;
Fs = 30000;
fftpts=1024*8;

[ERBforward,ERBfeedback,cf] = ERBfilt3(Fs, numChannels, lowFreq);
y = filtbank(ERBforward,ERBfeedback,noise);
subplot(2,1,1);
response = 20*log10(abs(fft(y',fftpts)));
freqScale = (0:fftpts-1)/fftpts*Fs;
axis([2 4 -70 10]);
semilogx(freqScale(1:fftpts/2),response(1:fftpts/2,:));
ylabel ('dB magnitude')
xlabel ('Frequency Hz')
title ('Noise passed through gammatone filter bank')


% ramp and scale the sounds
ramppts = rampdur * Fs;
stimpts = duration * Fs;
[rows,cols]=size(y);
ramp = maxval*[0:1/(ramppts -1):1 ones(1,(stimpts - 2*(ramppts ))) 1:-1/(ramppts -1):0];
ramp = ones(rows,1) * ramp;
for i = 1:rows
   m = abs(min(y(i,ramppts+1 : cols-ramppts-1 )));
   if (m < max(y(i,ramppts+1 : cols-ramppts-1 )))
      m = max(y(i, ramppts+1 : cols-ramppts-1 ));
   end
   y(i,:) = y(i,:)/m;
end

y2 = y .* ramp + 2048;

% do xcorrelations
% subplot(2,1,2)
% xcRef = max(xcorr(y2(14,:),y2(14,:)));
%for i = 1:20
%   temp(i) = max(xcorr(y2(i,:),y2(14,:)))/xcRef;
%end

%plot(cf,temp)
%ylabel ('relative correlation')
%xlabel ('Center Frequency Hz')
%title ('Correlation of Gammatone filtered noises')


% write to file
for i = [50:5:95]
   fname = ['GAM' num2str(cf(i)) fsuffix];
   fid = fopen([fpath fname],'wb');
   fwrite(fid,y2(i,:),'float');
   fclose(fid);
   disp(['filenum ' num2str(i) ' has max val : 'num2str(max(y2(i,:)))]);
   disp ([' and min value ' num2str(min(y2(i,:)))]);
end

   
   