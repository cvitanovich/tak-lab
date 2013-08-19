%Exit_Space3: Callback for Exit_ILDfreq pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig); H.rasterfig =[]; end
if(exist1('H.finalspikeratefig')) close(H.finalspikeratefig); H.finalspikeratefig =[]; end
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

set(H.buildOnly,'value',0)
XStimParams.buildOnly = 0;

%Save all values into current XStimParams
XStimParams.locations = [];
XStimParams.curr_ITD = str2num(get(H.space3_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.space3_ABL,'String'));
XStimParams.test_ISI = str2num(get(H.space3_ISI,'String'));
XStimParams.numreps = str2num(get(H.space3_numreps,'String'));

% deselect UseLastLocations_flag
GUI.UseLastLocations_flag =0;

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

close
%clear
clear s* t*
clear GLOBAL H

eval(['load ' FN.current_path 'H_xstim;'])
update_XStim
