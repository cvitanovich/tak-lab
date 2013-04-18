function zero_play_buffers
% flushes all play buffers with zeros
% TDT.bufpts -- no. buffer its for each buffer (1 x nBuffers vector)
% TDT.nPlayChannels -- no. of play channels
% TDT.nBuffers -- # of stimulus buffers in the play sequence
% TDt.stim_buffers -- buffer id values for each buffer

global TDT

S232('dropall');

for ch=1:TDT.nPlayChannels
	for buf=1:TDT.nBuffers
		S232('dpush',TDT.bufpts(buf));
		S232('value',0);
		S232('qpop16',TDT.stim_buffers(buf));
	end
end