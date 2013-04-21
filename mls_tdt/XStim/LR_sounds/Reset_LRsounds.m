%Reset_LRsounds: Callback for Reset_LRsounds pushbutton

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
   S232('trash');
   S232('dropall');
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig); H.rasterfig =[]; end
if(exist1('H.finalspikeratefig')) close(H.finalspikeratefig); H.finalspikeratefig =[]; end

set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

dur(1) = XStimParams.curr_stimdur;
dur(2) = XStimParams.curr_stimdur;

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
set(H.resetLRsounds,'Value',0);
XStimParams.reset_flag = 1;