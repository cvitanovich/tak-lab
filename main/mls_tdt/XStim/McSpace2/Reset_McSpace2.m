%reset McSpace2

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.McSpace2_finalspikeratefig')) close(H.McSpace2_finalspikeratefig);  H.McSpace2_finalspikeratefig =[]; end

GUI.locations1 = [];
setinfo_spacePicker

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

if(exist('H.McSpace2_finalspikeratefig')) close(H.McSpace2_finalspikeratefig); end

%Save all values into current XStimParams
XStimParams.locations = [];
if get(H.stim_type,'Value') ~= 9			% not from file
	XStimParams.bandwidth = str2num(get(H.McSpace2_bandwidth','String'));
	XStimParams.bandlimits(1) = str2num(get(H.McSpace2_bandlimits(1),'String'));
	XStimParams.bandlimits(2) = str2num(get(H.McSpace2_bandlimits(2),'String'));
end   
XStimParams.curr_ABL = str2num(get(H.McSpace2_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.McSpace2_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.McSpace2_ISI,'String'));
XStimParams.numreps = str2num(get(H.McSpace2_numreps,'String'));

% deselect UseLastLocations_flag
set(H.UseLastLocations,'value',0);
GUI.UseLastLocations_flag =0;
    
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetMcSpace2,'Value',0);
XStimParams.reset_flag = 1;

set(H.McSpace2_status,'BackgroundColor','blue');
set(H.McSpace2_status,'String','Status: reset, awaiting instructions');
