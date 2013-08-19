% SetInfo_Adapt2

% uses .eq HRIRs, therefore make sure earphone filters are also selected
   XStimParams.ephone_flag = 0;
   set(H.ephoneuseit,'Value',0);


    XStimParams.space_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
    XStimParams.ABLalone_flag = 0;
    
if get(H.Adapt2_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.Adapt2_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.Adapt2_DynamicRangepb,'BackgroundColor','red');
    set(H.Adapt2_DynamicRangeTxt,'visible','on');
else
    XStimParams.HiDynamicRange = 0;
    set(H.Adapt2_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.Adapt2_DynamicRangepb,'BackgroundColor','yellow');
    set(H.Adapt2_DynamicRangeTxt,'visible','off');
end

%%% values for SPLs
XStimParams.epoch_duration(1) = str2num(get(H.Adapt2_epoch_duration(1),'string'));
XStimParams.epoch_duration(2) = str2num(get(H.Adapt2_epoch_duration(2),'string'));

XStimParams.SPL1 = str2num(get(H.Adapt2_SPL1,'string'));
XStimParams.SPL2 = str2num(get(H.Adapt2_SPL2,'string'));

if min1(XStimParams.SPL2) < -10
    h = helpdlg('-10 dB is minimum allowed','SPL2 out of range');
    set(h,'position',[500 100 300 80]);
elseif min1(XStimParams.SPL2) == -10
    h = helpdlg('-10 dB is minimum allowed','SPL2 will be 0 dB');
    set(h,'position',[500 100 300 80]);
end
if max1(XStimParams.SPL2) > 31
    h = helpdlg('31 dB is maximum allowed','SPL2 out of range');
    set(h,'position',[500 100 300 80]);
end

if min1(XStimParams.SPL1) < -10
    h = helpdlg('-10 dB is minimum allowed','SPL1 out of range');
    set(h,'position',[500 100 300 80]);
elseif min1(XStimParams.SPL1) == -10
    h = helpdlg('-10 dB is minimum allowed','SPL1 will be 0 dB');
    set(h,'position',[500 100 300 80]);
end
if max1(XStimParams.SPL1) > 31
    h = helpdlg('31 dB is maximum allowed','SPL1 out of range');
    set(h,'position',[500 100 300 80]);
end

% update stimulus duration
XStimParams.silence_lead = str2num(get(H.Adapt2_silence_dur(1),'string'));
XStimParams.silence_mid = str2num(get(H.Adapt2_silence_dur(2),'string'));
XStimParams.silence_trail = str2num(get(H.Adapt2_silence_dur(3),'string'));

temp = length(XStimParams.SPL1);
XStimParams.curr_stimdur = temp*(length(XStimParams.SPL2) * sum(XStimParams.epoch_duration))+(temp-1)*XStimParams.silence_mid;
set(H.Adapt2_DUR, 'String',num2str(XStimParams.curr_stimdur));

if XStimParams.curr_stimdur > 36000
    h = helpdlg(['Total stimulus duration of ' num2str(XStimParams.curr_stimdur) ' excedes limit of 36000 in available memory'],'MEMORY LIMITATION EXCEDED');
    set(h,'position',[500 100 300 80]);
end

% update stimulus ABL
XStimParams.curr_ABL = str2num(get(H.Adapt2_ABL,'string'));

% update offsets
XStimParams.offset_el(2)= str2num(get(H.offset_el,'string'));
XStimParams.offset_az(2) = str2num(get(H.offset_az,'string'));

% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
    [FN.space_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space_std);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
end

if get(H.spaceflag,'value')
    XStimParams.space_flag = 1;
    set(H.spaceflag,'string','push for ABLequal FC')
    set(H.spaceflag,'BackgroundColor','yellow')
    set(H.spacefilepb,'string','Select FC Filter File');   
    set(H.spacefilehdr,'String', 'Fully-cued HRTF File:');
    set(H.spacefile,'string',FN.space_std);
    XStimParams.ABLequal_flag = 0;
    % check for FC file
    if get(H.spacefilepb,'Value')
        [FN.space_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',FN.space_std);
        end
        set(H.spacefilepb,'Value',0);
        FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
    end
    
else
    XStimParams.space_flag = 0;
    set(H.spaceflag,'string','push for Fully-cued')
    set(H.spaceflag,'BackgroundColor','red')
    set(H.spacefilepb,'String','Select ABLequal FC Filter File');   
    set(H.spacefilehdr,'String', 'ABLequal FC HRTF File:');
    set(H.spacefile,'string',FN.ablequal_std);
    XStimParams.ABLequal_flag = 1;
    
    % check for ABLequal file
    if get(H.spacefilepb,'Value')
        [FN.ablequal_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',FN.ablequal_std);
        end
        set(H.spacefilepb,'Value',0);
        FN.HRTFfiletype(7,2) = testHRTFfiletype(FN.space_path, FN.ablequal_std);
    end
end

% bandwidth & limits
XStimParams.bandlimits(1) = str2num(get(H.Adapt2_bandlimits1(1),'String'));
XStimParams.bandlimits(2) = str2num(get(H.Adapt2_bandlimits2(1),'String'));
XStimParams.bandlimits(3) = str2num(get(H.Adapt2_bandlimits1(2),'String'));
XStimParams.bandlimits(4) = str2num(get(H.Adapt2_bandlimits2(2),'String'));

% set numreps
XStimParams.numreps = str2num(get(H.Adapt2_numreps,'String'));

%Set space PickerFlag
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
if(exist1('H.Adapt2fig') & get(H.Adapt2_inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.Adapt2_inc_testnum,'Value',0);
end
set(H.Adapt2_data_FN,'String', FN.data);


% update stim_type1
stim_type = get(H.Adapt2_stim_type(1),'String');
stim_val = get(H.Adapt2_stim_type(1),'Value');

if ~strcmp(XStimParams.stim_type, deblank(stim_type(stim_val,:))) | ...
        get(H.Adapt2_stimFNpb(1),'Value')		% stim_type has changed or pb pushed
    XStimParams.stim_type = deblank(stim_type(stim_val,:));
    set(H.stim_type,'Value',stim_val);
    if(strcmp('File',deblank(XStimParams.stim_type)))
        [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
        set(H.Adapt2_stim_FN(1),'String', FN.stim);
        set(H.stim_filename,'String', FN.stim);
    end
    %setinfo_Xstim;
    set(H.Adapt2_stimFNpb(1),'Value',0);
end

%control visibility of FN vs freq
set(H.Adapt2_stimFNpb(1),'Visible','on');
set(H.Adapt2_stimFN(1),'Visible','on');
set(H.Adapt2_freqtxt(1),'Visible','off');
set(H.Adapt2_freq(1),'Visible','off');    
set(H.Adapt2_bandlimitstxt(1),'visible','off');
set(H.Adapt2_bandlimits1(1),'visible','off');
set(H.Adapt2_bandlimits2(1),'visible','off');

switch get(H.stim_type,'Value')
    case {1 2}
        set(H.Adapt2_freqtxt(1),'visible','on');
        set(H.Adapt2_freq(1),'visible','on');
        XStimParams.curr_freq = str2num(get(H.Adapt2_freq(1),'string'));
    case 6
        set(H.Adapt2_bandlimitstxt(1),'visible','on');
        set(H.Adapt2_bandlimits1(1),'visible','on');
        set(H.Adapt2_bandlimits2(1),'visible','on');
    case 8
    otherwise
    set(H.stim_type,'Value',8);
    disp('Unsupported stimulus type for Adapt2 test; reset to BROADBAND');
end

% update stim_type2
stim_type = get(H.Adapt2_stim_type(2),'String');
stim_val = get(H.Adapt2_stim_type(2),'Value');

if ~strcmp(XStimParams.stim_type2,deblank(stim_type(stim_val,:))) | ...
        get(H.Adapt2_stimFNpb(2),'Value')		% stim_type2 has changed
    XStimParams.stim_type2 = deblank(stim_type(stim_val,:));
    set(H.stim_type2,'Value',stim_val);
    if(strcmp('File',deblank(XStimParams.stim_type2)))
        [stim_from_file2, FN.stim2, FN.stim_path2] = Get_Stimulus_File(XStimParams, FN.stim_path2);
        set(H.Adapt2_stim_FN(2),'String', FN.stim2);
        set(H.stim_filename2,'String', FN.stim2);
    end
    set(H.Adapt2_stimFNpb(2),'Value',0);
end
%control visibility of FN vs freq
set(H.Adapt2_stimFNpb(2),'Visible','on');
set(H.Adapt2_stimFN(2),'Visible','on');
set(H.Adapt2_freqtxt(2),'Visible','off');
set(H.Adapt2_freq(2),'Visible','off');    
set(H.Adapt2_bandlimitstxt(2),'visible','off');
set(H.Adapt2_bandlimits1(2),'visible','off');
set(H.Adapt2_bandlimits2(2),'visible','off');

switch get(H.stim_type2,'Value')
    case {1 2}
        set(H.Adapt2_freqtxt(2),'visible','on');
        set(H.Adapt2_freq(2),'visible','on');
        XStimParams.curr_freq2 = str2num(get(H.Adapt2_freq(2),'string'));
    case 6
        set(H.Adapt2_bandlimitstxt(2),'visible','on');
        set(H.Adapt2_bandlimits1(2),'visible','on');
        set(H.Adapt2_bandlimits2(2),'visible','on');
    case 8
    otherwise
    set(H.stim_type2,'Value',8);
    disp('Unsupported stimulus type for Adapt2 test; reset to BROADBAND');
end





% check out modulation parameters for stim
XStimParams.mod_depth(1) = str2num(get(H.Adapt2_mod_depth,'String'));
mod_type = get(H.Adapt2_mod_type,'String');
mod_num = get(H.Adapt2_mod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Sq wave'
        set(H.Adapt2_mod_txt,'String', 'Freq:');
        set(H.Adapt2_mod_freq,'Visible','on');
        set(H.Adapt2_mod_phase,'Visible','on');
        set(H.Adapt2_mod_txtA,'Visible','on');
        set(H.Adapt2_mod_depth,'Visible','on');
        set(H.Adapt2_mod_txtB,'Visible','on');
        set(H.Adapt2_mod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Adapt2_mod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Adapt2_mod_phase,'String'));
    case 'Tone'
        set(H.Adapt2_mod_txt,'String', 'Freq:');
        set(H.Adapt2_mod_freq,'Visible','on');
        set(H.Adapt2_mod_phase,'Visible','on');
        set(H.Adapt2_mod_txtA,'Visible','on');
        set(H.Adapt2_mod_depth,'Visible','on');
        set(H.Adapt2_mod_txtB,'Visible','on');
        set(H.Adapt2_mod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Adapt2_mod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Adapt2_mod_phase,'String'));
    case 'File'
        set(H.Adapt2_mod_pb,'Visible','on');
        set(H.Adapt2_mod_freq,'Visible','off');
        set(H.Adapt2_mod_txtA,'Visible','on');
        set(H.Adapt2_mod_depth,'Visible','on');
        set(H.Adapt2_mod_phase,'Visible','off');
        set(H.Adapt2_mod_txtB,'Visible','off');
        if get(H.Adapt2_mod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.Adapt2_mod_pb,'Value',0);
        end
        set(H.Adapt2_mod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.Adapt2_mod_txt,'String', 'CutOff Freq:');
        set(H.Adapt2_mod_freq,'Visible','on');
        set(H.Adapt2_mod_txtA,'Visible','on');
        set(H.Adapt2_mod_depth,'Visible','on');
        set(H.Adapt2_mod_phase,'Visible','off');
        set(H.Adapt2_mod_pb,'Visible','off');
        set(H.Adapt2_mod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Adapt2_mod_freq,'String'));
        
    case 'None'
        set(H.Adapt2_mod_txt,'String', 'no mod  ');
        set(H.Adapt2_mod_freq,'Visible','off');
        set(H.Adapt2_mod_phase,'Visible','off');
        set(H.Adapt2_mod_pb,'Visible','off');
        set(H.Adapt2_mod_txtB,'Visible','off');
        set(H.Adapt2_mod_txtA,'Visible','off');
        set(H.Adapt2_mod_depth,'Visible','off');
        
    otherwise
end

update_XStim;

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])