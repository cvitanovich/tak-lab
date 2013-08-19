%SetABLStep

low_abl  = str2num(get(H.lowabl,'String'));
high_abl = str2num(get(H.highabl,'String'));
if(high_abl > -20)
   high_abl = -20;
   set(H.highabl,'String',num2str(high_abl));
end
numabls = str2num(get(H.numabls,'String'));

set(H.stepabl,'String',round(num2str((high_abl-low_abl)/(numabls-1))));