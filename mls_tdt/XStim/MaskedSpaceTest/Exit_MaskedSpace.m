%Exit_MaskedSpace: Callback for Exit_MaskedSpace pushbutton
reset_PD1;
m110dx( C_.STOP);

if(exist1('H.MaskedSpace_finalspikeratefig')) close(H.MaskedSpace_finalspikeratefig); H.MaskedSpace_finalspikeratefig =[]; end
if exist1('H.pickerfig') close(H.pickerfig); H.pickerfig = [];	end

set(H.locfile,'Enable','on');
set(H.locAZ,'Enable','on');
set(H.locEL,'Enable','on');
set(H.locuseit,'Enable','on');

%Save all values into current XStimParams
XStimParams.locations = [];
XStimParams.curr_stimdur = str2num(get(H.MaskedSpace_probeDUR,'String'));
XStimParams.curr_stimdur2 = str2num(get(H.MaskedSpace_maskDUR,'String'));
XStimParams.numreps = str2num(get(H.MaskedSpace_numreps,'String'));
XStimParams.test_ISI = 100;
set(H.Xstim_test_ISI,'string',num2str(XStimParams.test_ISI));

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
    if findstr(char(NAMES(i)),'maskedspace')
        ind = [ind i];
    end
    if findstr(char(NAMES(i)),'MaskedSpace')
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