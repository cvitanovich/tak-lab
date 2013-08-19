%Reset_space4

%if(S232('getAPlockstatus') & S232('getXBlockstatus'))
S232('PA4mute',1);
S232('PA4mute',2);
S232('PD1stop',1);
S232('PD1clear',1);
S232('trash');
S232('dropall');
%end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.space4_finalspikeratefig')) close(H.space4_finalspikeratefig);  H.space4_finalspikeratefig =[]; end
%if exist1('H.pickerfig') clearspace4axisCB;	end
GUI.locations1 = [];
setinfo_spacePicker
if exist1('H.space4_wts_fig') close(H.space4_wts_fig); H.space4_wts_fig = [];   end
if exist1('H.IAresponse_fig')   close(H.IAresponse_fig); H.IAresponse_fig=[];   end

set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

if(exist('H.space4_finalspikeratefig')) close(H.space4_finalspikeratefig); end

XStimParams.buildOnly = 0;

%Save all values into current XStimParams
XStimParams.locations = [];
XStimParams.curr_ITD = str2num(get(H.space4_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.space4_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.space4_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space4_ISI,'String'));
XStimParams.numreps = str2num(get(H.space4_numreps,'String'));

% deselect UseLastLocations_flag
set(H.UseLastLocations,'value',0);
GUI.UseLastLocations_flag =0;

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetspace4,'Value',0);
XStimParams.reset_flag = 1;

set(H.space4_status,'BackgroundColor','blue');
set(H.space4_status,'String','Status: reset, awaiting instructions');
