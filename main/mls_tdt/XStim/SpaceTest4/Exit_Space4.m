%Exit_space4: Callback for Exit_space4 pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
    S232('PA4mute',1);
    S232('PA4mute',2);
    S232('PD1stop',1);
    S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.space4_finalspikeratefig')) close(H.space4_finalspikeratefig); H.space4_finalspikeratefig =[]; end
if exist1('H.pickerfig') close(H.pickerfig); H.pickerfig = [];	end
if exist1('H.wts_fig') close(H.wts_fig);  H.wts_fig = [];  end
if exist1('H.IAresponse_fig')   close(H.IAresponse_fig); H.IAresponse_fig=[];   end

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

%Save all values into current XStimParams
XStimParams.HiDynamicRange = 0;
XStimParams.locations = [];
XStimParams.curr_ITD = str2num(get(H.space4_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.space4_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.space4_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space4_ISI,'String'));
XStimParams.numreps = str2num(get(H.space4_numreps,'String'));

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
    if findstr(char(NAMES(i)),'space4')
        ind = [ind i];
    end
    if findstr(char(NAMES(i)),'space4')
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