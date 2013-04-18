function TDT_attens(attens)
% for setting attenuations for the TDT
global TDT

for j=1:TDT.nPlayChannels
	S232('PA4atten',j,attens(j));
end