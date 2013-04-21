%Exit_Search: Callback for Exit_Search pushbutton

close (H.searchfig)
clear H.searchfig
clear H.psdaxes
if(exist1('H.spikefig'))
   close (H.spikefig)
   clear H.spikefig
   clear H.spikeaxes
end
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');

%clear
clear s* t*

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PD1stop',1);
   S232('PD1clear',1);
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('trash');
   S232('dropall');
   S232('S2close') % close application and release AP2 and XBUS locks
end

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

clear GLOBAL H
eval(['load ' FN.current_path 'H_xstim;'])