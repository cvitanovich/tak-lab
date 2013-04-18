function PD1_init
% generic function to prepare PD1 for single channel or multi-channel play
% TDT.din -- TDT device id no.
% TDT.Fs -- sampling rate in Hz
% TDT.npts_total_play -- total play its (all buffers)
% TDT.nPlayChannels -- number of channels to use
% TDT.nRecChannels -- no. of recording channels
global TDT

SRATE = 1e6 / TDT.Fs;
S232('PD1clear',TDT.din);
S232('PD1srate',TDT.din,SRATE);
S232('PD1npts',TDT.din,TDT.npts_total_play);
S232('PD1resetDSP',TDT.din,hex2dec('FFF'));
S232('dropall');
S232('PD1clrsched',TDT.din);
S232('PD1nstrms',TDT.din,TDT.nPlayChannels,TDT.nRecChannels);

% playback routing setup
for j=1:TDT.nPlayChannels
	CHAN = j - 1;
	S232('PD1addsimp', TDT.din, S232('DSPout', CHAN), S232('DAC', CHAN));
	S232('PD1specIB', TDT.din, S232('IB', CHAN), S232('DSPin', CHAN));
end

% record routing setup
if TDT.nRecChannels
	for k=1:TDT.nRecChannels
		CHAN = k - 1;
		S232('PD1specOB',TDT.din, S232('OB', CHAN), S232('ADC', CHAN));
	end
end
% LED thresholds
S232('PD1setIO',TDT.din,0.01,9.99,0.01,9.99);