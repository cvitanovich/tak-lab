%ILDAloneFileDiagBox
%Dialog Box to ask user for ildalone file

[FN.ildalone_eq, FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select VAS Filter File');

if(FN.ildalone ~= 0)
   set(H.ildalonefile,'String',FN.ildalone_eq);
	eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
end
