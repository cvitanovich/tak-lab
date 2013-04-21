%SetILDStep

low_ild  = str2num(get(H.lowild,'String'));
high_ild = str2num(get(H.highild,'String'));
numilds = str2num(get(H.numilds,'String'));

set(H.stepild,'String',num2str((high_ild-low_ild)/(numilds-1)));