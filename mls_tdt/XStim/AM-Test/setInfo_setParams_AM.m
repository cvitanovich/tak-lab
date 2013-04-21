function [] = setInfo_setParams_AM;

global GUI
global H
global FN
global XStimParams
global TDT

for count2 = 1:XStimParams.com_num_ted
	for count1 = 1:XStimParams.com_len_ted
		if XStimParams.amp_mat_ted(count2, 1, count1) ~= str2num(get(H.amp_mat_ted(count2, 1, count1),'String'))
            XStimParams.amp_mat_ted(count2, 1, count1) = str2num(get(H.amp_mat_ted(count2, 1, count1),'String'));
		end
		if XStimParams.amp_mat_ted(count2, 2, count1) ~= str2num(get(H.amp_mat_ted(count2, 2, count1),'String'))
            XStimParams.amp_mat_ted(count2, 2, count1) = str2num(get(H.amp_mat_ted(count2, 2, count1),'String'));
		end
        if XStimParams.amp_mat_ted(count2, 3, count1) ~= str2num(get(H.amp_mat_ted(count2, 3, count1),'String'))
            XStimParams.amp_mat_ted(count2, 3, count1) = str2num(get(H.amp_mat_ted(count2, 3, count1),'String'));
		end
        if XStimParams.amp_mat_ted(count2, 4, count1) ~= str2num(get(H.amp_mat_ted(count2, 4, count1),'String'))
            XStimParams.amp_mat_ted(count2, 4, count1) = str2num(get(H.amp_mat_ted(count2, 4, count1),'String'));
		end
        if XStimParams.len_mat_ted(count2, 1, count1) ~= str2num(get(H.len_mat_ted(count2, 1, count1),'String'))
            XStimParams.len_mat_ted(count2, 1, count1) = str2num(get(H.len_mat_ted(count2, 1, count1),'String'));
		end
        if XStimParams.len_mat_ted(count2, 2, count1) ~= str2num(get(H.len_mat_ted(count2, 2, count1),'String'))
            XStimParams.len_mat_ted(count2, 2, count1) = str2num(get(H.len_mat_ted(count2, 2, count1),'String'));
		end
        if XStimParams.len_mat_ted(count2, 3, count1) ~= str2num(get(H.len_mat_ted(count2, 3, count1),'String'))
            XStimParams.len_mat_ted(count2, 3, count1) = str2num(get(H.len_mat_ted(count2, 3, count1),'String'));
		end
        if XStimParams.len_mat_ted(count2, 4, count1) ~= str2num(get(H.len_mat_ted(count2, 4, count1),'String'))
            XStimParams.len_mat_ted(count2, 4, count1) = str2num(get(H.len_mat_ted(count2, 4, count1),'String'));
		end
	end
end

	
