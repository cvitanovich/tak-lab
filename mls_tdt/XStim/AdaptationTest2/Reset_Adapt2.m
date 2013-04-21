%reset Adapt2

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.Adapt2_finalspikeratefig')) close(H.Adapt2_finalspikeratefig);  H.Adapt2_finalspikeratefig =[]; end
if(exist1('H.Adapt2_epochFigure')) close(H.Adapt2_epochFigure);  H.Adapt2_epochFigure =[]; end

GUI.locations1 = [];
setinfo_spacePicker

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

if(exist('H.Adapt2_finalspikeratefig')) close(H.Adapt2_finalspikeratefig); end

%Save all values into current XStimParams
XStimParams.locations = [];
if get(H.stim_type,'Value') ~= 9			% not from file
	XStimParams.bandlimits(1) = str2num(get(H.Adapt2_bandlimits1(1),'String'));
	XStimParams.bandlimits(2) = str2num(get(H.Adapt2_bandlimits2(1),'String'));
end   
if get(H.stim_type2,'Value') ~= 9			% not from file
	XStimParams.bandlimits(3) = str2num(get(H.Adapt2_bandlimits1(2),'String'));
	XStimParams.bandlimits(4) = str2num(get(H.Adapt2_bandlimits2(2),'String'));
end   

XStimParams.curr_ABL = str2num(get(H.Adapt2_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.Adapt2_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.Adapt2_ISI,'String'));
XStimParams.numreps = str2num(get(H.Adapt2_numreps,'String'));

% deselect UseLastLocations_flag
set(H.UseLastLocations,'value',0);
GUI.UseLastLocations_flag =0;
    
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetAdapt2,'Value',0);
XStimParams.reset_flag = 1;

set(H.Adapt2_status,'BackgroundColor','blue');
set(H.Adapt2_status,'ForegroundColor','white');
set(H.Adapt2_status,'String','Status: reset, awaiting instructions');
