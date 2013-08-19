%hrtfFileDiagBox
%Dialog Box to ask user for earphone file

cd e:\spezio\matlab\scripts\tdt\XStim\HRTFfilts
[hrtffilename,hrtfpathname] = uigetfile('*.mat','Select HRTFFile for Analysis');

if(hrtfpathname ~= 0)
   set(hhrtffile,'String',[hrtfpathname hrtffilename]);
end

cd e:\spezio\matlab\scripts\tdt\XStim