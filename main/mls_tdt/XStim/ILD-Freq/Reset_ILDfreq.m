%Exit_ILDfreq: Callback for Exit_ILDfreq pushbutton

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
XStimParams.lofreq = str2num(get(H.lowfreq,'String'));
XStimParams.hifreq = str2num(get(H.highfreq,'String'));
XStimParams.numfreqs = str2num(get(H.numfreqs,'String'));
XStimParams.loild = str2num(get(H.lowild,'String'));
XStimParams.hiild = str2num(get(H.highild,'String'));
XStimParams.numilds = str2num(get(H.numilds,'String'));
XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetildfreq,'Value',0);
clear H.rasterfig