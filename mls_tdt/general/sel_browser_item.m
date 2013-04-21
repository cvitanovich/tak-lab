%Sel_Browser_Item
%Required for Browser usage
%Selects a browser item and makes it the string

browser_string = get(hbrowser,'String');
string_num = get(hbrowser,'Value');
chosen_string = browser_string(string_num,:);

if(isdir([pwd '\' chosen_string]))
   eval(['cd ' pwd '\' chosen_string])
end

