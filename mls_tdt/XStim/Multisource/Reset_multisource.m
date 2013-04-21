global GUI
global H
global XStimParams
global FN


%Exit_ILDfreq: Callback for Exit_ILDfreq pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist2('H.rasterfig'))
    try
        close(H.rasterfig);
    catch
    end
    H.rasterfig = [];
end
if(exist2('H.finalspikeratefig'))
    try
        close(H.finalspikeratefig);
    catch
    end
    H.finalspikeratefig = [];
end
if(exist2('H.pickerfig'))
    try
        close(H.pickerfig);
    catch
    end
    H.pickerfig = [];
end
if(exist2('H.setparamsfig'))
    try
        close(H.setparamsfig);
    catch
    end
    H.setparamsfig = [];
end

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
end

set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

% change status bar
set(H.status,'BackgroundColor','blue');
set(H.status,'String','Status: reset, awaiting instructions');

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetTwo_source,'Value',0);
XStimParams.reset_flag = 1;

