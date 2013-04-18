function [] = m_FF2

%function [] = m_FF2
% makes the frequency filter coefficients for
% narrow band freq tuning
% April 28, 2000
% may need to change wflag

wflag = 1;
minF = 2000;
maxF = 10000;
incF = 100;
Fs = 30000;
ntaps = 20;
outpath = 'd:\kip\coefs\';
outfname = 'filtcoef.003';

if wflag
   % set params for header
   n_lines = ntaps;
	nfreqs = (maxF - minF)/incF +1;;             % number of filters
	filetype = 2;                      				% = 1 for ints; 2 for floats
	n_channels = 2 * nfreqs;           				% for iir filtering (b & a coeffs)
	first_freq = minF;
	last_freq = maxF;
	comment1 = 'Gammatone filter coefficients for narrow band frequency tuning';
	comment2 = 'created with m_FF2.m, chan1 = b (forward) coeffs,  chan2 = a (feedback) coeffs';

% calculate info_blocks
	n_byte = nfreqs*2*4;
	info_blocks = ceil(n_byte/256)+2;

	mtlwh([outpath outfname], filetype, info_blocks, n_channels, n_lines, Fs,...
   	first_freq, last_freq, nfreqs, comment1, comment2);

 % write frequency matrix
	dblfreq = [minF:incF:maxF; minF:incF:maxF];   % just to match HRTF dir matrix
	mtlwdir([outpath outfname], dblfreq);
end

i = 0;
for freq = minF:incF:maxF
   i = i+1;
   tempfor = zeros(1,ntaps);
   tempback = zeros(1,ntaps);
   
   [forward, feedback] = make1_ERBfilt(freq,Fs,0);
   tempfor(1:length(forward)) = forward;
   tempback(1:length(feedback)) = feedback;
   
	if wflag
	  status = mtlwch([outpath outfname],tempfor, i*2-1);
	  status = mtlwch([outpath outfname],tempback, i*2);
	end
end

fclose('all');
