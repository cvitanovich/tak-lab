%Exit_space_RIR: Callback for Exit_space_RIR pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
    S232('PA4mute',1);
    S232('PA4mute',2);
    S232('PD1stop',1);
    S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.finalspikerasterfig')) close(H.finalspikerasterfig);  H.finalspikerasterfig =[]; end

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
    S232('trash');
    S232('dropall');
    S232('S2close') % close application and release AP2 and XBUS locks
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');
set(H.locfile,'Enable','on');
set(H.locAZ,'Enable','on');
set(H.locEL,'Enable','on');
set(H.locuseit,'Enable','on');

XStimParams.ephone_flag = 0;
set(H.ephoneuseit,'Value',0);

XStimParams.buildOnly = 0;

%Save all values into current XStimParams
XStimParams.HiDynamicRange = 0;
XStimParams.locations = [];
XStimParams.curr_stimdur = str2num(get(H.space_RIR_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space_RIR_ISI,'String'));
XStimParams.numreps = str2num(get(H.space_RIR_numreps,'String'));

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

% deselect UseLastLocations_flag
GUI.UseLastLocations_flag =0;

close
clear s* t* GLOBAL

ind = [];
NAMES = fieldnames(H);
for i = 1:length(NAMES)
    if findstr(char(NAMES(i)),'space_RIR')
        ind = [ind i];
    end
    if findstr(char(NAMES(i)),'space_rir')
        ind = [ind i];
    end
    if findstr(char(NAMES(i)),'picker')
        ind = [ind i];
    end
end

for i = 1:length(ind)
    if isfield(H,NAMES(ind(i)))
        H = rmfield(H,NAMES(ind(i)));
    end
end
clear NAMES ind i

XStimParams.reset_flag = 1;
if exist1('H.composite_fig')
    setinfo_composite;
end
update_XStim