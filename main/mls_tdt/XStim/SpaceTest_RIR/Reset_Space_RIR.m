%Reset_space_RIR

%if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
   S232('trash');
   S232('dropall');
%end
m110dx( C_.STOP);

if(exist1('H.finalspikerasterfig')) close(H.finalspikerasterfig);  H.finalspikerasterfig =[]; end

set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

%Save all values into current XStimParams
XStimParams.locations = [];
XStimParams.curr_stimdur = str2num(get(H.space_RIR_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space_RIR_ISI,'String'));
XStimParams.numreps = str2num(get(H.space_RIR_numreps,'String'));
    
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetspace_RIR,'Value',0);
XStimParams.reset_flag = 1;

set(H.space_RIR_status,'BackgroundColor','blue');
set(H.space_RIR_status,'String','Status: reset, awaiting instructions');
