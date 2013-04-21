%Reset_Composite

% may not have established process, because, unlike the other modules,
% composite does not set up a process when the GUI is made.
try(S232('getS2primary'));
catch
    if(S232('S2init', 0, 'INIT_PRIMARY', 20000) == 0)
        disp('Cannot initialize a primary process')
        return;
    end
end

% get locks
if(S232('APlock', 100, 0) == 0)
    disp('Cannot acquire lock on AP2 Card')
    s232('S2close');
    return;
end
if(S232('XBlock', 100, 0) == 0)
    disp('Cannot acquire lock on X Bus')
    s232('APunlock', 0);
    s232('S2close');
    return;
end

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
    S232('PA4mute',1);
    S232('PA4mute',2);
    S232('PD1stop',1);
    S232('PD1clear',1);
    
end
m110dx( C_.STOP);

% eliminate earlier tests
if isfield(H,'composite_test_type_txt')    
    if Test_num
        for itest = 1: Test_num 
            delete(H.composite_test_type_txt(itest));
            delete(H.composite_dataFN(itest));
        end
    end
    H.composite_test_type_txt = [];
    H.composite_dataFN = [];
end
Test_num = 0;
F.str_loc = F.str_locStart;

if exist1('H.composite_datafig')
    for ifig = 1:length( H.composite_datafig)
        close(H.composite_datafig(ifig));
    end
end
if isfield(H,'composite_datafig')   
    H.composite_datafig = [];
    H.composite_spikeaxes = [];
end

% reset duration
XStimParams.curr_stimdur = 100;
set(H.composite_DUR, 'String', num2str(100))

set(H.composite_num_stims_txt,'string',['no stimuli chosen']);

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
    S232('trash');
    S232('dropall');
    S232('S2close') % close application and release AP2 and XBUS locks
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

%Save all values into current XStimParams
XStimParams.curr_ABL = str2num(get(H.composite_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.composite_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.composite_ISI,'String'));
XStimParams.numreps = str2num(get(H.composite_numreps,'String'));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.composite_reset,'Value',0);
XStimParams.reset_flag = 0;

clear Test_array
XStimParams.reset_flag = 0;
set(H.composite_status,'BackgroundColor','blue');
set(H.composite_status,'String','Status: reset, awaiting instructions');