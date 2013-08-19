%SetBandWidth

XStimParams.bandwidth(1) = str2num(get(H.bandwidth,'String'));
if exist1(H.bandwidth2)
   XStimParams.bandwidth(2) = str2num(get(H.bandwidth2,'String'));
else
   XStimParams.bandwidth(2) = 1/3;
end
%Make sure earphone filters are also selected
XStimParams.ephone_flag = 1;
set(H.ephoneuseit,'Value',1);
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
return
