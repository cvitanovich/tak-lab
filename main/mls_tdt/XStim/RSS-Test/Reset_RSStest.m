%Exit_ABLtest: Callback for Exit_ABLtest pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist('H.rasterfig')) close(H.rasterfig); end
if(exist('H.finalspikeratefig')) close(H.finalspikeratefig); end

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

%Save all values into current XStimParams
XStimParams.loabl = str2num(get(H.lowabl,'String'));
XStimParams.hiabl = str2num(get(H.highabl,'String'));
XStimParams.numabls = str2num(get(H.numabls,'String'));
XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;']);
eval(['save ' FN.current_path 'FN_current FN;']);
set(H.resetabltest,'Value',0);
clear H.rasterfig