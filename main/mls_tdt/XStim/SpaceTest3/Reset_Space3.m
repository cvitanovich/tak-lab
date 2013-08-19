%Reset_space3

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.finalspikeratefig')) close(H.finalspikeratefig);  H.finalspikeratefig =[]; end
%if exist1('H.pickerfig') clearspaceaxisCB;	end
GUI.locations1 = [];
setinfo_spacePicker
if exist1('H.wts_fig') close(H.wts_fig); H.wts_fig = [];   end
if exist1('H.IAresponse_fig')   close(H.IAresponse_fig); H.IAresponse_fig=[];   end

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

if(exist('H.finalspikeratefig')) close(H.finalspikeratefig); end

%Save all values into current XStimParams
XStimParams.locations = [];
XStimParams.curr_freq = str2num(get(H.space3_freq,'String'));
XStimParams.curr_ITD = str2num(get(H.space3_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.space3_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.space3_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space3_ISI,'String'));
XStimParams.numreps = str2num(get(H.space3_numreps,'String'));

% deselect UseLastLocations_flag
set(H.UseLastLocations,'value',0);
GUI.UseLastLocations_flag =0;
    
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetspace,'Value',0);
XStimParams.reset_flag = 1;

set(H.status,'BackgroundColor','blue');
set(H.status,'String','Status: reset, awaiting instructions');
