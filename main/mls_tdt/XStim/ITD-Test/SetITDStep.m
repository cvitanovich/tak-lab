%SetITDStep

low_itd  = str2num(get(H.lowitd,'String'));
high_itd = str2num(get(H.highitd,'String'));
if(high_itd > 250)
   high_itd = 250;
   set(H.highitd,'String',num2str(high_itd));
end
if(low_itd < -250)
   low_itd = -250;
   set(H.lowitd,'String',num2str(low_itd));
end

numitds = str2num(get(H.numitds,'String'));

set(H.stepitd,'String',num2str(round((high_itd-low_itd)/(numitds-1))));