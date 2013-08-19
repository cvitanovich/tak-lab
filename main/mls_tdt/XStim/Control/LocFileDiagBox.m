%LOCFileDiagBox
%Dialog Box to ask user for file for 1 specified location

[FN.loc, FN.loc_path] = uigetfile([FN.loc_path '*.*'],'Select Location Filter File');

if(FN.loc_path ~= 0)
    set(H.locfile,'String',[FN.loc_path FN.loc]);
	eval(['save ' FN.current_path 'XStimParams_current XStimParams']);
    eval(['save ', FN.current_path  'FN_current FN'])
end
