%Exit_ILDfreq: Callback for Exit_ILDfreq pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.space2_finalspikeratefig')) close(H.space2_finalspikeratefig);  H.space2_finalspikeratefig =[]; end
%if exist1('H.pickerfig') clearspaceaxisCB;	end
GUI.locations1 = [];
setinfo_spacePicker
if exist1('H.space2_wts_fig') close(H.space2_wts_fig); H.space2_wts_fig = [];   end
if exist1('H.IAresponse_fig')   close(H.IAresponse_fig); H.IAresponse_fig=[];   end

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

if(exist('H.space2_finalspikeratefig')) close(H.space2_finalspikeratefig); end

%set(H.space2_buildOnly,'value',0)
%XStimParams.buildOnly = 0;

%Save all values into current XStimParams
XStimParams.locations = [];
if get(H.stim_type,'Value') ~= 9			% not from file
	XStimParams.bandwidth = str2num(get(H.space2_bandwidth','String'));
	XStimParams.bandlimits(1) = str2num(get(H.space2_bandlimits1,'String'));
	XStimParams.bandlimits(2) = str2num(get(H.space2_bandlimits2,'String'));
end   
XStimParams.curr_ITD = str2num(get(H.space2_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.space2_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.space2_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space2_ISI,'String'));
XStimParams.numreps = str2num(get(H.space2_numreps,'String'));

% deselect UseLastLocations_flag
set(H.UseLastLocations,'value',0);
GUI.UseLastLocations_flag =0;
    
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetspace,'Value',0);
XStimParams.reset_flag = 1;

set(H.space2_status,'BackgroundColor','blue');
set(H.space2_status,'String','Status: reset, awaiting instructions');
