%Form_Directory_List
%Required for Browser usage
%Changes 'string' property to the current directory listing

dlist = dir;
set(hbrowser,'String',str2mat(dlist(1:length(dlist)).name));