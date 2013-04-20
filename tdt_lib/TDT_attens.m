function TDT_attens(TDT)
% for setting attenuations for the TDT
% TDT.attens = attenuation settings for each channel

for j=1:TDT.nPlayChannels
	S232('PA4atten',j,TDT.attens(j));
end