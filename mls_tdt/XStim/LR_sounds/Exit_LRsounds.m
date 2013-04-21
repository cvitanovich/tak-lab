%Exit_LRsounds

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
   S232('trash');
   S232('dropall');
   S232('S2close') % close application and release AP2 and XBUS locks
end
m110dx( C_.STOP);

if(exist1('H.rasterfig')) close(H.rasterfig); H.rasterfig =[]; end
if(exist1('H.finalspikeratefig')) close(H.finalspikeratefig); H.finalspikeratefig =[]; end
if(exist1('H.LRsoundsfig')) close(H.LRsoundsfig); H.LRsoundsfig =[]; end


set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');
set(H.locfile,'Enable','on');
set(H.locAZ,'Enable','on');
set(H.locEL,'Enable','on');
set(H.locuseit,'Enable','on');


XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
clear s* t* BUF itd* dur

clear GLOBAL H
eval(['load ' FN.current_path 'H_xstim;'])
temp = fieldnames(H);
for i = 1:length(temp)
    if findstr(temp{i},'LR')
        H = rmfield(H,temp{i});
    end
end
clear temp i
XStimParams.reset_flag = 1;