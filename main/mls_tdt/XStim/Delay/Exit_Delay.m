%Exit_2source: Callback for Exit pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig);  H.rasterfig =[]; end
if(exist1('H.finalspikeratefig')) close(H.finalspikeratefig); H.finalspikeratefig =[]; end
if exist1('H.pickerfig') close(H.pickerfig); H.pickerfig = [];	end

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
   S232('S2close')              % close application and release AP2 and XBUS locks
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');
set(H.locfile,'Enable','on');
set(H.locAZ,'Enable','on');
set(H.locEL,'Enable','on');
set(H.locuseit,'Enable','on');

%Save all values into current XStimParams
XStimParams.locations = [];

% XStimParams.bandwidth(1) = str2num(get(H.bandwidth1,'String'));
% XStimParams.bandlimits(1) = str2num(get(H.bandlimits1,'String'));
% XStimParams.bandlimits(2) = str2num(get(H.bandlimits2,'String'));
% %XStimParams.bandwidth(2) = str2num(get(H.bandwidth2,'String'));

XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

close
clear s* t*
clear GLOBAL H
eval(['load ' FN.current_path 'H_xstim;'])

if(findobj('Tag','Stimuli'))
    close('Stimuli');
end
if(findobj('Tag','Raster'))
    close('Raster');
end
if(findobj('Tag','Delay plot'))
    close('Delay plot');
end
if(findobj('Tag','Test Raster'))
    close('Test Raster');
end
