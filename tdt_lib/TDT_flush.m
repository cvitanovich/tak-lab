function TDT_flush(TDT)
% commands to stopTDT and flush all it's settings

for ch=1:TDT.nPlayChannels
	S232('PA4mute',ch);
end

S232('PD1stop',TDT.din);
S232('PD1clrIO',TDT.din);
S232('PD1clear',TDT.din);

S232('trash');
S232('dropall');

S232('APunlock',0);
S232('XBunlock',0);
S232('S2close');