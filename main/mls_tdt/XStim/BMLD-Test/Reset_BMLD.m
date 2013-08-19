%Exit_BMLD: Callback for Exit_BMLD pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.BMLDspikefig')) close(H.BMLDspikefig);  H.BMLDspikefig =[]; end

GUI.locations1 = [];

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

% change status bar
set(H.status,'BackgroundColor','blue');
set(H.status,'String','Status: reset, awaiting instructions');

%Save all values into current XStimParams
if get(H.stim_type,'Value') ~= 9			% not from file
	XStimParams.bandwidth = [str2num(get(H.bandwidth1,'String')) str2num(get(H.bandwidth2,'String'))];
	XStimParams.bandlimits(1) = str2num(get(H.bandlimits1,'String'));
	XStimParams.bandlimits(2) = str2num(get(H.bandlimits2,'String'));
end   
XStimParams.curr_ITD = str2num(get(H.ITD(1),'String'));
XStimParams.ITD2 = str2num(get(H.ITD(2),'String'));
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetBMLD,'Value',0);
XStimParams.reset_flag = 1;