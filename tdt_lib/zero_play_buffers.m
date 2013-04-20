function zero_play_buffers(TDT)
% flushes all play buffers with zeros
% TDT.bufpts -- no. buffer its for each buffer (1 x nBuffers vector)
% TDT.nPlayChannels -- no. of play channels
% TDT.nBuffers -- # of stimulus buffers in the play sequence
% TDT.stim_buffers -- buffer id values for each buffer

S232('dropall');

for ch=1:TDT.nPlayChannels
    nbuffers=length(TDT.playpts{ch});
	for buf=1:nbuffers
		S232('dpush',TDT.playpts{ch}(buf));
		S232('value',0);
		S232('qpop16',TDT.stim_buffers{ch}(buf));
	end
end