%setInfo_freq


% increment test number
if(exist1('H.freqtestfig') & get(H.inc_testnum,'Value'))
	XStimParams.testnum = XStimParams.testnum +1;
	set(H.testnum, 'String', num2str(XStimParams.testnum))
   update_dataFN;
   set(H.inc_testnum,'Value',0);
end
set(H.recorddata_FN,'String', FN.data);

low_freq  = str2num(get(H.lowfreq,'String'));
high_freq = str2num(get(H.highfreq,'String'));
numfreqs = str2num(get(H.numfreqs,'String'));

set(H.stepfreq,'String',num2str(round((high_freq-low_freq)/(numfreqs-1))));

stim_type = get(H.freqstim_type,'String');
stim_val = get(H.freqstim_type,'Value');
XStimParams.stim_type = deblank(stim_type(stim_val,:));
set(H.stim_type,'Value', stim_val);
clear stim_type

% save current parameter values
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim