%Reset_PlayFile

%if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
   S232('trash');
   S232('dropall');
%end
m110dx( C_.STOP);


if(exist('H.PlayFile_finalspikeratefig')) close(H.PlayFile_finalspikeratefig); end

%Save all values into current XStimParams
XStimParams.locations = [];
if get(H.stim_type,'Value') ~= 9			% not from file
	XStimParams.bandwidth = str2num(get(H.PlayFile_bandwidth','String'));
	XStimParams.bandlimits(1) = str2num(get(H.PlayFile_bandlimits1,'String'));
	XStimParams.bandlimits(2) = str2num(get(H.PlayFile_bandlimits2,'String'));
end   
XStimParams.curr_ITD = str2num(get(H.PlayFile_ITD,'String'));
XStimParams.curr_ILD = str2num(get(H.PlayFile_ILD,'String'));
XStimParams.curr_ABL = str2num(get(H.PlayFile_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.PlayFile_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.PlayFile_ISI,'String'));
XStimParams.numreps = str2num(get(H.PlayFile_numreps,'String'));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetPlayFile,'Value',0);
XStimParams.reset_flag = 1;

set(H.PlayFile_status,'BackgroundColor','blue');
set(H.PlayFile_status,'String','Status: reset, awaiting instructions');
