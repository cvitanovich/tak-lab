%dataFileDiagBox
%Dialog Box to ask user for earphone file

cd e:\spezio\data
[datafilename,datapathname] = uigetfile('*.mat','Select DataFile for Analysis');

if(datapathname ~= 0)
   set(hdatafile,'String',[datapathname datafilename]);
end
cd e:\spezio\matlab\scripts\tdt\XStim