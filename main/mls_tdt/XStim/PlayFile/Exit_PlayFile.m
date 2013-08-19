%Exit_PlayFile: Callback for Exit_PlayFile pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
    S232('PA4mute',1);
    S232('PA4mute',2);
    S232('PD1stop',1);
    S232('PD1clear',1);
    S232('trash');
    S232('dropall');
    S232('S2close') % close application and release AP2 and XBUS locks
end
m110dx( C_.STOP);

if(exist1('H.PlayFile_finalspikeratefig')) close(H.PlayFile_finalspikeratefig); H.PlayFile_finalspikeratefig =[]; end

%Save all values into current XStimParams
XStimParams.HiDynamicRange = 0;
XStimParams.locations = [];
XStimParams.curr_ITD = str2num(get(H.PlayFile_ITD,'String'));
XStimParams.curr_ILD = str2num(get(H.PlayFile_ILD,'String'));
XStimParams.curr_ABL = str2num(get(H.PlayFile_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.PlayFile_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.PlayFile_ISI,'String'));
XStimParams.numreps = str2num(get(H.PlayFile_numreps,'String'));

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

close
clear s* t* GLOBAL

ind = [];
NAMES = fieldnames(H);
for i = 1:length(NAMES)
    if findstr(char(NAMES(i)),'PlayFile')
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