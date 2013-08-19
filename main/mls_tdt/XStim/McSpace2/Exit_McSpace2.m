% Exit_McSpace2

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

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.McSpace2_finalspikeratefig')) close(H.McSpace2_finalspikeratefig); H.McSpace2_finalspikeratefig =[]; end
if exist1('H.pickerfig') close(H.pickerfig); H.pickerfig = [];	end

set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');
set(H.locfile,'Enable','on');
set(H.locAZ,'Enable','on');
set(H.locEL,'Enable','on');
set(H.locuseit,'Enable','on');

XStimParams.ephone_flag = 0;
set(H.ephoneuseit,'Value',0);

% if stim_type is BBN-LNN, switch to BBN
stim_type = get(H.McSpace2_stim_type,'String');
stim_val = get(H.McSpace2_stim_type,'Value');
if stim_val == 10   stim_val = 8;   end
set(H.stim_type,'Value',stim_val);
XStimParams.stim_type = deblank(stim_type(stim_val,:));

   %Save all values into current XStimParams
XStimParams.locations = [];
XStimParams.curr_stimdur = str2num(get(H.McSpace2_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.McSpace2_ISI,'String'));
XStimParams.numreps = str2num(get(H.McSpace2_numreps,'String'));
XStimParams.rndSPLs = 1;

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

XStimParams.reset_flag = 1;

% deselect UseLastLocations_flag
GUI.UseLastLocations_flag =0;

close
clear s* t* GLOBAL

ind = [];
NAMES = fieldnames(H);
for i = 1:length(NAMES)
    if findstr(char(NAMES(i)),'McSpace2') | findstr(char(NAMES(i)),'McSpace2') | findstr(char(NAMES(i)),'picker')
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