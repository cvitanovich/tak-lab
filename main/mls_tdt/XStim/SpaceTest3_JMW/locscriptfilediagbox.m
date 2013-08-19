%LocScriptFileDiagBox
%Dialog Box to ask user for Location Script file

[FN.script, FN.script_path] = uigetfile('*.scr','Select Location Script File');

if(FN.script ~= 0)
   set(H.locscriptfile,'String',[FN.script_path FN.scriptfile]);
	eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
end
