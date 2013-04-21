%Exit_ITDtest: Callback for Exit_ITDtest pushbutton

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

set(H.itd_buildOnly,'value',0)
XStimParams.buildOnly = 0;

%Save all values into current XStimParams
XStimParams.loitd = str2num(get(H.lowitd,'String'));
XStimParams.hiitd = str2num(get(H.highitd,'String'));
XStimParams.numitds = str2num(get(H.numitds,'String'));
XStimParams.curr_ABL = str2num(get(H.itd_ABL,'String'));
XStimParams.curr_ILD = str2num(get(H.itd_ILD,'String'));
XStimParams.curr_stimdur = str2num(get(H.itd_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.itd_ISI,'String'));
XStimParams.numreps = str2num(get(H.itd_numreps,'String'));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetitdtest,'Value',0);
clear H.rasterfig
XStimParams.reset_flag = 1;