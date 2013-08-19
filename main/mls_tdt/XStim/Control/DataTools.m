%DataTools: Callback for DataTools popupmenu

val = get(H.datatools,'Value');

switch val
   
case 1
    prep4_MLSweights;
case 2
   IAMultRegress;
otherwise
   disp('Not implemented')
   set(H.datatools,'Value',1);
end

return