function [] = m_FF
%function [] = m_FF
% makes the frequency filter coefficients for
% narrow band freq tuning
% January 28, 1998
% may need to change wflag

wflag = 0;

minF = 2000;
maxF = 10000;
Fs = 30000;
fftpts=1024*8;

outpath = 'd:\kip\temp\';
outfname = 'filtcoef.00x';

% first make 120 gammatone filters with 1/30 octave (see ERBfilt3.m) spacing from 1000 Hz upwards
[ERBforward,ERBfeedback,cf] = ERBfilt3(Fs,0,1000);

% plot the response to a square pulse
impulse = [1 zeros(1,1023)];
y = filtbank(ERBforward,ERBfeedback,impulse);
figure
subplot(2,1,1)
FT = abs(fft(y',fftpts));
response = 20*log10(FT);
freqScale = (0:fftpts-1)/fftpts*Fs;
axis([2 4 -70 10]);
semilogx(freqScale(1:fftpts/2),response(1:fftpts/2,:));
ylabel ('dB magnitude')
xlabel ('Frequency Hz')
title ('ERB gammatone filter bank')
sum(FT(274:3277,:))';   % just to look see

% plot to compare Q10dB pts
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


% pull out center frequencies of interest
[m,n] = size(response)
index = '';
for i = 1:n
   index(i) = find(response(1:fftpts/2,i) == max(response(1:fftpts/2,i)));
end
freq = freqScale(index)'                       % picks freqs at max response
index = find((freq >= minF) & (freq <= maxF)); % b/wn min and max
freq = round(freq(index));                     % 
if wflag

response = response(1:fftpts/2,index);        % 
ERBforward = ERBforward(index,:);
ERBfeedback = ERBfeedback(index,:);

% make sizes of each channel the same
[mfor,nfor] = size(ERBforward);
[mfeed,nfeed] = size(ERBfeedback);
if nfor > nfeed
   ERBfeedback = [ERBfeedback zeros(mfor, nfor-nfeed)];
   n_lines = nfor;           % # a or b coefficients per filter
elseif nfor < nfeed
   ERBforward = [ERBforward zeros(mfor, nfeed-nfor)];
   n_lines = nfeed;          % # a or b coefficients per filter
end


% set params for header
nfreqs = mfor;                     % number of filters
filetype = 2;                      % = 1 for ints; 2 for floats
n_channels = 2 * nfreqs;           % for iir filtering (b & a coeffs)
first_freq = min(freq);            %
last_freq = max(freq);
comment1 = 'Gammatone filter coefficients for narrow band frequency tuning';
comment2 = 'created with m_FF.m, chan1 = b coeffs,  chan2 = a coeffs';

% calculate info_blocks
n_byte = nfreqs*2*4;
info_blocks = ceil(n_byte/256)+2;

mtlwh([outpath outfname], filetype, info_blocks, n_channels, n_lines, Fs,...
   first_freq, last_freq, nfreqs, comment1, comment2);

 % write frequency matrix
dblfreq = [freq; freq];   % just to match HRTF dir matrix
mtlwdir([outpath outfname], dblfreq);

% write coefficients
for channel = 1:nfreqs
   display(['writing frequency filter number ' num2str(channel)]);
   % b coefficients:
  status = mtlwch([outpath outfname],ERBforward(channel,:), channel*2-1);
   % a coefficients:
  status = mtlwch([outpath outfname],ERBfeedback(channel,:), channel*2);
end

fclose('all');
end