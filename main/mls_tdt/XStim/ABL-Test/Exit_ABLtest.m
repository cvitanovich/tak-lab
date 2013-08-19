%Exit_ILDfreq: Callback for Exit_ILDfreq pushbutton

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
if(exist1('H.pickerfig'))   close(H.pickerfig); H.pickerfig = [];   end

set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');
set(H.locfile,'Enable','on');
set(H.locAZ,'Enable','on');
set(H.locEL,'Enable','on');
set(H.locuseit,'Enable','on');

set(H.abl_buildOnly,'value',0)
XStimParams.buildOnly = 0;

%Save all values into current XStimParams
XStimParams.loabl = str2num(get(H.lowabl,'String'));
XStimParams.hiabl = str2num(get(H.highabl,'String'));
XStimParams.numabls = str2num(get(H.numabls,'String'));
XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
XStimParams.curr_ILD = str2num(get(H.ILD,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.loc_flag = 0;

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

close
%clear
%clear s* t* BUF abl*
%clear GLOBAL H

ind = [];
NAMES = fieldnames(H);
for i = 1:length(NAMES)
    if findstr(char(NAMES(i)),'abl') | findstr(char(NAMES(i)),'ABL')
        ind = [ind i];
    end
end

for i = 1:length(ind)
    H = rmfield(H,NAMES(ind(i)));
end
clear NAMES ind i

%eval(['load ' FN.current_path 'H_xstim;'])

XStimParams.reset_flag = 1;
if exist1('H.composite_fig')
    setinfo_composite;
end
update_XStim
