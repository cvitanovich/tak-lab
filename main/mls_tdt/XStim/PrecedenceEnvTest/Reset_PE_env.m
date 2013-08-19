%Reset_PE_env

%if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
   %end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.finalspikerasterfig')) close(H.finalspikerasterfig);  H.finalspikerasterfig =[]; end
if exist1('H.PE_env_fig') close(H.PE_env_fig); H.PE_env_fig=[]; end

%if exist1('H.pickerfig') clearspaceaxisCB;	end
GUI.locations1 = [];
setinfo_spacePicker;

%if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
   %end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

if(exist('H.finalspikeratefig')) close(H.finalspikeratefig); end

set(H.recorddata_FN,'String',FN.data);
% change status bar
set(H.status,'BackgroundColor','blue');
set(H.status,'String','Status: reset, awaiting instructions');

%Save all values into current XStimParams
XStimParams.locations = [];
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
%XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));

% deselect UseLastLocations_flag
set(H.UseLastLocations,'value',0);
GUI.UseLastLocations_flag =0;

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
XStimParams.reset_flag = 1;