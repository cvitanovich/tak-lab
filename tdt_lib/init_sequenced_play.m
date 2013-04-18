function init_sequenced_play
% initialize sequenced play for the TDT
global TDT

S232('seqplay',TDT.play_spec);
if TDT.nRecChannels
	S232('seqrecord',TDT.rec_spec);
end
S232('PD1arm',TDT.din);
S232('pfireall');
S232('PD1go',TDT.din);