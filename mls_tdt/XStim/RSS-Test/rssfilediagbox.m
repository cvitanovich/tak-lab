%RSSFileDiagBox
%Dialog Box to ask user for RSS file

cd e:\spezio\stimuli

[FN.rss,FN.rss_path] = uigetfile('*.mat','Select Earphone Filter File');

if(rsspathname ~= 0)
   set(H.rssfile,'String',[FN.rss_path FN.rss]);
   eval(['save ' FN.current_path 'XStimParams_current XStimParams;']);
end

eval(['cd ' FN.home]);