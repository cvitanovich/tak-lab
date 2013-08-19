%Reset_MaskedSpace
reset_PD1;
m110dx( C_.STOP );

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.MaskedSpace_finalspikeratefig')) close(H.MaskedSpace_finalspikeratefig);  H.MaskedSpace_finalspikeratefig =[]; end
%if exist1('H.pickerfig') clearspaceaxisCB;	end
GUI.locations1 = [];
setinfo_spacePicker
set(H.pickerfig,'Name','Probe Space Picker')

%Save all values into current XStimParams
XStimParams.locations = [];
XStimParams.curr_stimdur = str2num(get(H.MaskedSpace_probeDUR,'String'));
XStimParams.curr_stimdur2 = str2num(get(H.MaskedSpace_maskDUR,'String'));
XStimParams.numreps = str2num(get(H.MaskedSpace_numreps,'String'));

% deselect UseLastLocations_flag
set(H.UseLastLocations,'value',0);
GUI.UseLastLocations_flag =0;

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetMaskedSpace,'Value',0);
XStimParams.reset_flag = 1;

set(H.MaskedSpace_status,'BackgroundColor','blue');
set(H.MaskedSpace_status,'String','Status: reset, awaiting instructions');
