%Reset_LongSAMs

m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.LongSAMs_finalspikeratefig')) close(H.LongSAMs_finalspikeratefig);  H.LongSAMs_finalspikeratefig =[]; end
if(exist1('H.LongSAMs_PERIODfig')) close(H.LongSAMs_PERIODfig);  H.LongSAMs_PERIODfig =[]; end
%if exist1('H.pickerfig') clearspaceaxisCB;	end
GUI.locations1 = [];
setinfo_spacePicker

set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

if(exist('H.LongSAMs_finalspikeratefig')) close(H.LongSAMs_finalspikeratefig); end

XStimParams.buildOnly = 0;

%Save all values into current XStimParams
XStimParams.locations = [];
XStimParams.curr_ABL = str2num(get(H.LongSAMs_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.LongSAMs_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.LongSAMs_ISI,'String'));
XStimParams.numreps = str2num(get(H.LongSAMs_numreps,'String'));

% deselect UseLastLocations_flag
set(H.UseLastLocations,'value',0);
GUI.UseLastLocations_flag =0;
    
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetLongSAMs,'Value',0);
XStimParams.reset_flag = 1;

set(H.LongSAMs_status,'BackgroundColor','blue');
set(H.LongSAMs_status,'String','Status: reset, awaiting instructions');
