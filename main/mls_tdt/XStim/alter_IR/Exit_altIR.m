%Exit_altIR: Callback for Exit_altIR pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.altIR_finalspikeratefig')) close(H.altIR_finalspikeratefig); H.altIR_finalspikeratefig =[]; end
if exist1('H.pickerfig') close(H.pickerfig); H.pickerfig = [];	end

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

set(H.altIR_buildOnly,'value',0)
XStimParams.buildOnly = 0;

%Save all values into current XStimParams
XStimParams.locations = [];
XStimParams.curr_ITD = str2num(get(H.altIR_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.altIR_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.altIR_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.altIR_ISI,'String'));
XStimParams.numreps = str2num(get(H.altIR_numreps,'String'));

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

XStimParams.reset_flag = 1;

% deselect UseLastLocations_flag
GUI.UseLastLocations_flag =0;

close
clear s* t*

ind = [];
NAMES = fieldnames(H);
for i = 1:length(NAMES)
    if findstr(char(NAMES(i)),'altIR')
        ind = [ind i];
    end
    if findstr(char(NAMES(i)),'alter')
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