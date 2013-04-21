function [] = setInfo_setParams;

% SetInfo_spacePicker
global GUI
global H
global FN
global XStimParams
global TDT

for count1 = 1:size(XStimParams.sub_locations, 1)
	% update 
	
	XStimParams.stim_type_ted(count1) = get(H.Twosrc_stim_type(count1), 'Value');
    
	if XStimParams.stim_type_ted(count1) ~= H.stim_type(count1) | get(H.Twosrc_stim_filenamepb(count1),'Value')
		H.stim_type(count1) = XStimParams.stim_type_ted(count1);
		if H.stim_type(count1) == 9;
			[stim_from_file, tempstr1, FN.stim_path] = Get_Stimulus_File_ted(char(XStimParams.stimfile(count1)), FN.stim_path);
			XStimParams.stimfile(count1) = cellstr(tempstr1);
            FN.stim_ted(count1) = cellstr(tempstr1);
            set(H.Twosrc_stim_filename(count1), 'String', char(FN.stim_ted(count1)));
		end
		%setinfo_Xstim;
		set(H.Twosrc_stim_filenamepb(count1),'Value',0);
	end

	if XStimParams.freq(count1) ~= str2num(get(H.Twosrc_freq(count1),'String'))
		XStimParams.freq(count1) = str2num(get(H.Twosrc_freq(count1),'String'));
	end
	
	if XStimParams.factor(count1) ~= str2num(get(H.factor(count1),'String'))
        XStimParams.factor(count1) = str2num(get(H.factor(count1),'String'));
	end
	
	if XStimParams.bandwidth(count1) ~= str2num(get(H.bandwidth(count1),'String'))
        XStimParams.bandwidth(count1) = str2num(get(H.bandwidth(count1),'String'));
    end
    
    if XStimParams.bandlimits1(count1) ~= str2num(get(H.bandlimits1(count1),'String'))
		XStimParams.bandlimits1(count1) = str2num(get(H.bandlimits1(count1),'String'));
    end
	
	if XStimParams.bandlimits2(count1) ~= str2num(get(H.bandlimits2(count1),'String'))
		XStimParams.bandlimits2(count1) = str2num(get(H.bandlimits2(count1),'String'));
	end
  
    %control visibility of stim params
	switch XStimParams.stim_type_ted(count1)
		case {1, 2}
			set(H.Twosrc_stim_filenamepb(count1), 'Enable','off');
			set(H.Twosrc_stim_filename(count1), 'String','');
% 			set(H.nostim_txt(count1),'Visible','off');
			set(H.freq_txt(count1),'Visible','on');
			set(H.Twosrc_freq(count1), 'Enable','on');       
			set(H.factor_txt(count1),'Visible','on');
			set(H.factor(count1), 'Enable','on');
			set(H.bandwidth_txt(count1),'Visible','on');
			set(H.bandwidth(count1), 'Enable','off');
			set(H.limit_txt(count1),'Visible','on');
			set(H.bandlimits1(count1), 'Enable','off');
			set(H.bandlimits2(count1), 'Enable','off');
		case 9
			set(H.Twosrc_stim_filenamepb(count1), 'Enable','on');
			set(H.Twosrc_stim_filename(count1), 'String', char(XStimParams.stimfile(count1)));
% 			set(H.nostim_txt(count1),'Visible','off');
			set(H.freq_txt(count1),'Visible','on');
			set(H.Twosrc_freq(count1), 'Enable','off');       
			set(H.factor_txt(count1),'Visible','on');
			set(H.factor(count1), 'Enable','off');
			set(H.bandwidth_txt(count1),'Visible','on');
			set(H.bandwidth(count1), 'Enable','off');
			set(H.limit_txt(count1),'Visible','on');
			set(H.bandlimits1(count1), 'Enable','off');
			set(H.bandlimits2(count1), 'Enable','off');
		case 10
			set(H.Twosrc_stim_filenamepb(count1), 'Enable','off');
			set(H.Twosrc_stim_filename(count1), 'String','');
% 			set(H.nostim_txt(count1),'Visible','off');
			set(H.freq_txt(count1),'Visible','on');
			set(H.Twosrc_freq(count1), 'Enable','off');       
			set(H.factor_txt(count1),'Visible','on');
			set(H.factor(count1), 'Enable','off');
			set(H.bandwidth_txt(count1),'Visible','on');
			set(H.bandwidth(count1), 'Enable','off');
			set(H.limit_txt(count1),'Visible','on');
			set(H.bandlimits1(count1), 'Enable','off');
			set(H.bandlimits2(count1), 'Enable','off');
		otherwise 
			set(H.Twosrc_stim_filenamepb(count1), 'Enable','off');
			set(H.Twosrc_stim_filename(count1), 'String','');
% 			set(H.nostim_txt(count1),'Visible','off');
			set(H.freq_txt(count1),'Visible','on');
			set(H.Twosrc_freq(count1), 'Enable','on');       
			set(H.factor_txt(count1),'Visible','on');
			set(H.factor(count1), 'Enable','on');
			set(H.bandwidth_txt(count1),'Visible','on');
			set(H.bandwidth(count1), 'Enable','on');
			set(H.limit_txt(count1),'Visible','on');
			set(H.bandlimits1(count1), 'Enable','on');
			set(H.bandlimits2(count1), 'Enable','on');
	end
    
% 	XStimParams.mod_type(count1) = get(H.Twosrc_mod_type(count1), 'Value');
% 	
% 	
% 	if get(H.Twosrc_mod_pb(count1),'Value')
% 		[mod_from_file(count1), FN.mod(count1), FN.mod_path(count1)] = Get_Stimulus_File(XStimParams, FN.mod_path(count1));
% 		set(H.Twosrc_mod_pb(count1),'Value',0);
% 	end
% 	
	
	XStimParams.mod_type_ted(count1) = get(H.Twosrc_mod_type(count1), 'Value');

	if XStimParams.mod_type_ted(count1) ~= H.mod_type(count1) | get(H.Twosrc_mod_pb(count1),'Value')
		H.mod_type(count1) = XStimParams.mod_type_ted(count1);
		if H.mod_type(count1) == 3;
		    [mod_from_file, tempstr1, FN.mod_path] = Get_Stimulus_File_ted(char(XStimParams.modfile(count1)), FN.mod_path);
			XStimParams.modfile(count1) = cellstr(tempstr1);
            FN.mod_ted(count1) = cellstr(tempstr1);
            set(H.mod_file_txt(count1), 'String', char(FN.mod_ted(count1)));
		end
		set(H.Twosrc_mod_pb(count1),'Value',0);
	end
	
% 	if H.mod_num(count1) ~= XStimParams.mod_type(count1)
% 		H.mod_num(count1) = XStimParams.mod_type(count1);
% 	end
	
	if XStimParams.mod_freq(count1) ~= str2num(get(H.Twosrc_mod_freq(count1),'String'));
		XStimParams.mod_freq(count1) = str2num(get(H.Twosrc_mod_freq(count1),'String'));
	end
	
	if XStimParams.mod_depth(count1) ~= str2num(get(H.Twosrc_mod_depth(count1),'String'));
		XStimParams.mod_depth(count1) = str2num(get(H.Twosrc_mod_depth(count1),'String'));
	end
	
	if XStimParams.mod_phase(count1) ~= str2num(get(H.Twosrc_mod_phase(count1),'String'));
		XStimParams.mod_phase(count1) = str2num(get(H.Twosrc_mod_phase(count1),'String'));
	end
		
	%control visibility of mod params
	switch XStimParams.mod_type_ted(count1)
		case 1
			set(H.Twosrc_mod_pb(count1),'Enable','off');
			set(H.mod_file_txt(count1), 'String', '');
% 			H.nomod_txt(count1)
			set(H.mod_freq_txt(count1), 'Visible', 'on');
			set(H.Twosrc_mod_freq(count1),'Enable','on');
			set(H.mod_depth_txt, 'Visible', 'on');
			set(H.Twosrc_mod_depth(count1),'Enable','on');
			set(H.mod_phase_txt, 'Visible', 'on');
			set(H.Twosrc_mod_phase(count1),'Enable','on');
		case 2
			set(H.Twosrc_mod_pb(count1),'Enable','off');
			set(H.mod_file_txt(count1), 'String', '');
% 			H.nomod_txt(count1)
			set(H.mod_freq_txt(count1), 'Visible', 'on');
			set(H.Twosrc_mod_freq(count1),'Enable','off');
			set(H.mod_depth_txt, 'Visible', 'on');
			set(H.Twosrc_mod_depth(count1),'Enable','on');
			set(H.mod_phase_txt, 'Visible', 'on');
			set(H.Twosrc_mod_phase(count1),'Enable','off');

		case 3
			set(H.Twosrc_mod_pb(count1),'Enable','on');
			set(H.mod_file_txt(count1), 'String', char(XStimParams.modfile(count1)));
% 			H.nomod_txt(count1)
			set(H.mod_freq_txt(count1), 'Visible', 'on');
			set(H.Twosrc_mod_freq(count1),'Enable','off');
			set(H.mod_depth_txt, 'Visible', 'on');
			set(H.Twosrc_mod_depth(count1),'Enable','off');
			set(H.mod_phase_txt, 'Visible', 'on');
			set(H.Twosrc_mod_phase(count1),'Enable','off');
		case 4
			set(H.Twosrc_mod_pb(count1),'Enable','off');
			set(H.mod_file_txt(count1), 'String', '');
% 			H.nomod_txt(count1)
			set(H.mod_freq_txt(count1), 'Visible', 'on');
			set(H.Twosrc_mod_freq(count1),'Enable','off');
			set(H.mod_depth_txt, 'Visible', 'on');
			set(H.Twosrc_mod_depth(count1),'Enable','off');
			set(H.mod_phase_txt, 'Visible', 'on');
			set(H.Twosrc_mod_phase(count1),'Enable','off');
		otherwise
	end
end
 
if XStimParams.curr_stimdur ~= str2num(get(H.DUR,'String'))
    XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
end

isthereafile_flag = 0;
for count1 = 1:size(XStimParams.sub_locations, 1)
    if H.stim_type(count1) == 9
        set(H.DUR,'Enable','off');
        clear dir
        D = dir([FN.stim_path char(XStimParams.stimfile(count1))]);
        XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
        set(H.DUR,'string',num2str(XStimParams.curr_stimdur));
        isthereafile_flag = 1;    
    end
end
if ~isthereafile_flag
    set(H.DUR,'Enable','on');
end


