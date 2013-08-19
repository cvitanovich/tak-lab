%SetBandLimits

XStimParams.bandlimits(1) = str2num(get(H.bandlimits1,'String'));
XStimParams.bandlimits(2) = str2num(get(H.bandlimits2,'String'));
%Make sure earphone filters are also selected
XStimParams.ephone_flag = 1;
set(H.ephoneuseit,'Value',1);
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
return
