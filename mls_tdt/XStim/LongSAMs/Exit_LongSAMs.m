%Exit_LongSAMs: Callback for Exit_LongSAMs pushbutton

m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.LongSAMs_finalspikeratefig')) close(H.LongSAMs_finalspikeratefig); H.LongSAMs_finalspikeratefig =[]; end
if(exist1('H.LongSAMs_PERIODfig')) close(H.LongSAMs_PERIODfig);  H.LongSAMs_PERIODfig =[]; end
if exist1('H.pickerfig') close(H.pickerfig); H.pickerfig = [];	end

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
XStimParams.curr_ABL = str2num(get(H.LongSAMs_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.LongSAMs_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.LongSAMs_ISI,'String'));
XStimParams.numreps = str2num(get(H.LongSAMs_numreps,'String'));

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);
XStimParams.mod_type = 'None';
XStimParams.curr_stimdur = 100;

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

% deselect UseLastLocations_flag
GUI.UseLastLocations_flag =0;

close
clear s* t* GLOBAL

ind = [];
NAMES = fieldnames(H);
for i = 1:length(NAMES)
    if findstr(char(NAMES(i)),'LongSAMs')
        ind = [ind i];
    end
    if findstr(char(NAMES(i)),'LongSAMs')
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