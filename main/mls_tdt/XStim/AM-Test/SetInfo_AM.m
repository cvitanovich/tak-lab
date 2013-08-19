% SetInfo_AM
% combines all the other 'set...' scripts into one


if get(H.spaceflag,'Value') & (XStimParams.space_flag == 0)
    XStimParams.space_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
elseif get(H.ildaloneflag,'Value') & (XStimParams.ildalone_flag == 0)
    XStimParams.space_flag = 0;
    XStimParams.ildalone_flag = 1;
    XStimParams.itdalone_flag = 0;
elseif get(H.itdaloneflag,'Value') & (XStimParams.itdalone_flag == 0)
    XStimParams.space_flag = 0;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 1;
end


set(H.spaceflag,'Value',XStimParams.space_flag);
set(H.ildaloneflag,'value',XStimParams.ildalone_flag);
set(H.itdaloneflag,'value',XStimParams.itdalone_flag);

% control visibility of file picking for ILA, ITA or fully-cued
if get(H.spaceflag,'Value')
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');   
    set(H.ildalonefilehdr,'Enable','off');
    set(H.ildalonefile,'Enable','off');
    set(H.ildalonepb,'Enable','off');  
    set(H.itdalonefilehdr,'Enable','off');
    set(H.itdalonefile,'Enable','off');
    set(H.itdalonepb,'Enable','off');  
    set(H.ITDchk,'Value',0);
    set(H.ILDchk,'Value',0);
elseif get(H.ildaloneflag,'Value')
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');
    set(H.ildalonefilehdr,'Enable','on');
    set(H.ildalonefile,'Enable','on');
    set(H.ildalonepb,'Enable','on');
    set(H.itdalonefilehdr,'Enable','off');
    set(H.itdalonefile,'Enable','off');
    set(H.itdalonepb,'Enable','off');  
    set(H.ITDchk,'Value',1);
    set(H.ILDchk,'Value',0);
elseif get(H.itdaloneflag,'Value')
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');
    set(H.ildalonefilehdr,'Enable','off');
    set(H.ildalonefile,'Enable','off');
    set(H.ildalonepb,'Enable','off');
    set(H.itdalonefilehdr,'Enable','on');
    set(H.itdalonefile,'Enable','on');
    set(H.itdalonepb,'Enable','on');  
    set(H.ITDchk,'Value',0);
    set(H.ILDchk,'Value',1);
end

% check for ILA pushbutton (file-picker)
if get(H.ildalonepb,'Value')
    [FN.ildalone,FN.ILA_path] = uigetfile([FN.ILA_path '*.flt'],'Select ILD-alone HRTF File');
    if(FN.ILA_path ~= 0)
        set(H.ildalonefile,'String',FN.ildalone);
    end
    set(H.ildalonepb,'Value',0);
    FN.HRTFfiletype(2) = testHRTFfiletype(FN.ILA_path, FN.ildalone);
end

% check for ITA pushbutton (file-picker)
if get(H.itdalonepb,'Value')
    [FN.itdalone,FN.ITA_path] = uigetfile([FN.ITA_path '*.ita'],'Select ITD-alone HRTF File');
    if(FN.ITA_path ~= 0)
        set(H.itdalonefile,'String',FN.itdalone);
    end
    set(H.itdalonepb,'Value',0);
    FN.HRTFfiletype(3) = testHRTFfiletype(FN.ITA_path, FN.itdalone);
end

% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
    [FN.space,FN.space_path] = uigetfile([FN.space_path '*.eq'],'Select Fully-cued HRTF File');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1) = testHRTFfiletype(FN.space_path, FN.space);
end

%SetPickerFlag
XStimParams.picker_flag = 1;
if exist1('H.pickerfig') 	
    if isempty(H.pickerfig)
        XStimParams.locations = [];
        H.pickerfig = spacePicker;
    end
else
    XStimParams.locations = [];
    H.pickerfig = spacePicker;  
end

% increment test number
if(exist1('H.AM_fig') & get(H.inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.inc_testnum,'Value',0);
end

% update stim control area
if get(H.space3_stim_type,'Value') == 9								% from file
    if get(H.space3_stimFNpb,'Value')
        [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
        set(H.space3_stimFNpb,'Value',0);
	end
    set(H.space3_stimFN,'String',FN.stim);
    set(H.stim_filename,'String',FN.stim);

    set(H.stim_filename,'Enable','on');
    set(H.space3_stimFN,'Enable','on');
    set(H.space3_stimFNpb,'Enable','on');
    set(H.stim_filename,'Visible','on');
    set(H.space3_stimFN,'Visible','on');
    set(H.space3_stimFNpb,'Visible','on');
   
    set(H.com_len_ted,'Enable','off');
    set(H.com_len_ted_txt,'Enable','off');
    set(H.space3_stim_seg_pb,'Enable','off');
    set(H.com_num_ted,'Enable','off');
    set(H.com_len_ted_txt2,'Enable','off');
    set(H.com_len_ted,'Visible','off');
    set(H.com_len_ted_txt,'Visible','off');
    set(H.space3_stim_seg_pb,'Visible','off');    
    set(H.com_num_ted,'Visible','off');
    set(H.com_len_ted_txt2,'Visible','off');
    
    D = dir([FN.stim_path FN.stim]);
    XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.DUR, 'String', num2str(XStimParams.curr_stimdur));
elseif get(H.space3_stim_type,'Value') == 11	
    set(H.stim_filename,'Enable','off');
    set(H.space3_stimFN,'Enable','off');
    set(H.space3_stimFNpb,'Enable','off');
    set(H.stim_filename,'Visible','off');
    set(H.space3_stimFN,'Visible','off');
    set(H.space3_stimFNpb,'Visible','off');
   
    set(H.com_len_ted,'Enable','on');
    set(H.com_len_ted_txt,'Enable','on');
    set(H.space3_stim_seg_pb,'Enable','on');
    set(H.com_num_ted,'Enable','on');
    set(H.com_len_ted_txt2,'Enable','on');
    set(H.com_len_ted,'Visible','on');
    set(H.com_len_ted_txt,'Visible','on');
    set(H.space3_stim_seg_pb,'Visible','on');    
    set(H.com_num_ted,'Visible','on');
    set(H.com_len_ted_txt2,'Visible','on');
    
    count1 = str2num(get(H.com_len_ted, 'String'));
    count2 = str2num(get(H.com_num_ted, 'String'));
    if H.com_len_ted > XStimParams.com_len_ted
        XStimParams.amp_mat_ted(1:XStimParams.com_num_ted, :, (XStimParams.com_len_ted + 1):count1) = ones(XStimParams.com_num_ted, 4, count1 - XStimParams.com_len_ted);
        XStimParams.len_mat_ted(1:XStimParams.com_num_ted, :, (XStimParams.com_len_ted + 1):count1) = ones(XStimParams.com_num_ted, 4, count1 - XStimParams.com_len_ted);
    end
    if H.com_num_ted > XStimParams.com_num_ted
        XStimParams.amp_mat_ted((XStimParams.com_num_ted + 1):count2, :, :) = ones(count2 - XStimParams.com_num_ted, 4, count1);
        XStimParams.len_mat_ted((XStimParams.com_num_ted + 1):count2, :, :) = ones(count2 - XStimParams.com_num_ted, 4, count1);
    end
    XStimParams.com_len_ted = str2num(get(H.com_len_ted, 'String')); 
    XStimParams.com_num_ted = str2num(get(H.com_num_ted, 'String')); 
    if get(H.space3_stim_seg_pb,'Value')
        setParams_AM;
    end
else  
    set(H.stim_filename,'Enable','off');
    set(H.space3_stimFN,'Enable','off');
    set(H.space3_stimFNpb,'Enable','off');
    set(H.stim_filename,'Visible','off');
    set(H.space3_stimFN,'Visible','off');
    set(H.space3_stimFNpb,'Visible','off');
    
    set(H.com_len_ted,'Enable','off');
    set(H.com_len_ted_txt,'Enable','off');
    set(H.space3_stim_seg_pb,'Enable','off');
    set(H.com_num_ted,'Enable','off');
    set(H.com_len_ted_txt2,'Enable','off');
    set(H.com_len_ted,'Visible','off');
    set(H.com_len_ted_txt,'Visible','off');
    set(H.space3_stim_seg_pb,'Visible','off');    
    set(H.com_num_ted,'Visible','off');
    set(H.com_len_ted_txt2,'Visible','off');
    
end




%update display of data filename
set(H.space3_recorddata_FN,'String', FN.data);

% check out modulation parameters for stim
XStimParams.mod_depth(1) = str2num(get(H.Spacetest3mod_depth,'String'));
mod_type = get(H.Spacetest3mod_type,'String');
mod_num = get(H.Spacetest3mod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));
set(H.Spacetest3mod_txtC,'Visible','off');
set(H.Spacetest3mod_freq,'Visible','off');
set(H.Spacetest3mod_phase,'Visible','off');
set(H.Spacetest3mod_pb,'Visible','off');
set(H.Spacetest3mod_txtB,'Visible','off');
set(H.Spacetest3mod_txtA,'Visible','off');
set(H.Spacetest3mod_depth,'Visible','off');

switch XStimParams.mod_type
    case 'Sq wave'
        set(H.Spacetest3mod_txt,'String', 'Freq:');
        set(H.Spacetest3mod_freq,'Visible','on');
        set(H.Spacetest3mod_phase,'Visible','on');
        set(H.Spacetest3mod_txtA,'Visible','on');
        set(H.Spacetest3mod_depth,'Visible','on');
        set(H.Spacetest3mod_txtB,'Visible','on');
        XStimParams.mod_freq(1) = str2num(get(H.Spacetest3mod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Spacetest3mod_phase,'String'));
        set(H.Spacetest3mod_txtC,'Visible','on');
        set(H.Spacetest3mod_txtC,'String','pi/2 phase begins at onset of squarewave max');
    case 'Tone'
        set(H.Spacetest3mod_txt,'String', 'Freq:');
        set(H.Spacetest3mod_freq,'Visible','on');
        set(H.Spacetest3mod_phase,'Visible','on');
        set(H.Spacetest3mod_txtA,'Visible','on');
        set(H.Spacetest3mod_depth,'Visible','on');
        set(H.Spacetest3mod_txtB,'Visible','on');
        XStimParams.mod_freq(1) = str2num(get(H.Spacetest3mod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Spacetest3mod_phase,'String'));
        set(H.Spacetest3mod_txtC,'Visible','on');
        set(H.Spacetest3mod_txtC,'String','pi phase begins at peak of sinewave');
    case 'File'
        set(H.Spacetest3mod_pb,'Visible','on');
        set(H.Spacetest3mod_txtA,'Visible','on');
        set(H.Spacetest3mod_depth,'Visible','on');
        if get(H.Spacetest3mod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.Spacetest3mod_pb,'Value',0);
        end
        set(H.Spacetest3mod_txt,'String', FN.mod);
    case 'LP Noise'
        set(H.Spacetest3mod_txt,'String', 'CutOff Freq:');
        set(H.Spacetest3mod_freq,'Visible','on');
        set(H.Spacetest3mod_txtA,'Visible','on');
        set(H.Spacetest3mod_depth,'Visible','on');
        XStimParams.mod_freq(1) = str2num(get(H.Spacetest3mod_freq,'String'));
    case 'None'
        set(H.Spacetest3mod_txt,'String', 'no mod  ');
    otherwise
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim