%SetFreqStep

low_freq  = str2num(get(H.lowfreq,'String'));
high_freq = str2num(get(H.highfreq,'String'));
numfreqs = str2num(get(H.numfreqs,'String'));

set(H.stepfreq,'String',num2str(round((high_freq-low_freq)/(numfreqs-1))));