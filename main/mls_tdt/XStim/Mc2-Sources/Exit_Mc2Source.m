%Exit_Mc2source: Callback for Exit pushbutton

%if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
   S232('trash');
   S232('dropall');
   S232('S2close')              % close application and release AP2 and XBUS locks
   %end
m110dx( C_.STOP);

if exist1('H.rasterfig') close(H.rasterfig);  H.rasterfig =[]; end
if exist1('H.finalspikeratefig') close(H.finalspikeratefig); H.finalspikeratefig =[]; end
if exist1('H.pickerfig') close(H.pickerfig); H.pickerfig = [];	end
if exist1('H.Mc2SourceRLFfig') close(H.Mc2SourceRLFfig); H.Mc2SourceRLFfig = []; end
if exist1('FigParams.handle') close(FigParams.handle); end
clear FigParams

set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');
set(H.locfile,'Enable','on');
set(H.locAZ,'Enable','on');
set(H.locEL,'Enable','on');
set(H.locuseit,'Enable','on');

%Save all values into current XStimParams
XStimParams.curr_stimdur = 100;
XStimParams.test_ISI = 100;
XStimParams.focalProb =.8;
XStimParams.HiDynamicRange = 0;
XStimParams.locations = [];
XStimParams.bandlimits(1) = str2num(get(H.bandlimits1,'String'));
XStimParams.bandlimits(2) = str2num(get(H.bandlimits2,'String'));
XStimParams.bandlimits(3) = str2num(get(H.bandlimits3,'String'));
XStimParams.bandlimits(4) = str2num(get(H.bandlimits4,'String'));

XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

close
clear s* t*
clear GLOBAL H
eval(['load ' FN.current_path 'H_xstim;'])
update_XStim