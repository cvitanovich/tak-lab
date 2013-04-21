%SpaceFileDiagBox3
%Dialog Box to ask user for space file

[FN.space,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select VAS Filter File');

if(FN.space ~= 0)
   set(H.spacefile,'String',[FN.space_path FN.space]);
   eval(['save ' FN.current_path 'FN_current FN;'])
end
