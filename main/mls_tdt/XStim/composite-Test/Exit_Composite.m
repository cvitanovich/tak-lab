%Exit_Composite: Callback for Exit_Composite pushbutton

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

% clear data figures
if exist1('H.composite_datafig')
    for ifig = 1:length( H.composite_datafig)
        close(H.composite_datafig(ifig));
    end
end
if isfield(H,'composite_datafig')   rmfield(H,'composite_datafig'); end

% reset duration
XStimParams.curr_stimdur = 100;

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
XStimParams.locations = [];
XStimParams.curr_ABL = str2num(get(H.composite_ABL,'String'));
XStimParams.curr_stimdur = 100;
XStimParams.test_ISI = str2num(get(H.composite_ISI,'String'));
XStimParams.numreps = str2num(get(H.composite_numreps,'String'));

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

close
%clear
clear s* t* GLOBAL H Test_array

eval(['load ' FN.current_path 'H_xstim;'])