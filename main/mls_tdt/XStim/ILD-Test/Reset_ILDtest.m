%Exit_ILDtest: Callback for Exit_ILDtest pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig); H.rasterfig =[]; end
if(exist1('H.finalspikeratefig')) close(H.finalspikeratefig); H.finalspikeratefig =[]; end

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

set(H.ild_buildOnly,'value',0)
XStimParams.buildOnly = 0;

%Save all values into current XStimParams
XStimParams.loild = str2num(get(H.lowild,'String'));
XStimParams.hiild = str2num(get(H.highild,'String'));
XStimParams.numilds = str2num(get(H.numilds,'String'));
XStimParams.curr_freq = str2num(get(H.ild_freq,'String'));
XStimParams.curr_ABL = str2num(get(H.ild_ABL,'String'));
XStimParams.curr_ITD = str2num(get(H.ild_ITD,'String'));
XStimParams.curr_stimdur = str2num(get(H.ild_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ild_ISI,'String'));
XStimParams.numreps = str2num(get(H.ild_numreps,'String'));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetildtest,'Value',0);
clear H.rasterfig
XStimParams.reset_flag = 1;