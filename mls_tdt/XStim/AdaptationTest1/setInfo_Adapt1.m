% SetInfo_Adapt1

% uses .eq HRIRs, therefore make sure earphone filters are NOT also selected
   XStimParams.ephone_flag = 0;
   set(H.ephoneuseit,'Value',0);


    XStimParams.space_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
    XStimParams.ABLalone_flag = 0;
    
if get(H.Adapt1_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.Adapt1_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.Adapt1_DynamicRangepb,'BackgroundColor','red');
    set(H.Adapt1_DynamicRangeTxt,'visible','on');
else
    XStimParams.HiDynamicRange = 0;
    set(H.Adapt1_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.Adapt1_DynamicRangepb,'BackgroundColor','yellow');
    set(H.Adapt1_DynamicRangeTxt,'visible','off');
end

%%% values for SPLs
XStimParams.epoch_duration(1) = str2num(get(H.Adapt1_epoch_duration(1),'string'));
XStimParams.epoch_duration(2) = str2num(get(H.Adapt1_epoch_duration(2),'string'));

XStimParams.SPL1 = str2num(get(H.Adapt1_SPL1,'string'));
XStimParams.SPL2 = str2num(get(H.Adapt1_SPL2,'string'));

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
XStimParams.silence_lead = str2num(get(H.Adapt1_silence_dur(1),'string'));
XStimParams.silence_mid = str2num(get(H.Adapt1_silence_dur(2),'string'));
XStimParams.silence_trail = str2num(get(H.Adapt1_silence_dur(3),'string'));

temp = length(XStimParams.SPL1);
XStimParams.curr_stimdur = temp*(length(XStimParams.SPL2) * sum(XStimParams.epoch_duration))+(temp-1)*XStimParams.silence_mid;
set(H.Adapt1_DUR, 'String',num2str(XStimParams.curr_stimdur));

if XStimParams.curr_stimdur > 36000
    h = helpdlg(['Total stimulus duration of ' num2str(XStimParams.curr_stimdur) ' excedes limit of 36000 in available memory'],'MEMORY LIMITATION EXCEDED');
    set(h,'position',[500 100 300 80]);
end

% update stimulus ABL
XStimParams.curr_ABL = str2num(get(H.Adapt1_ABL,'string'));

% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
    [FN.space_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space_std);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
end


% bandwidth & limits
XStimParams.bandlimits(1) = str2num(get(H.Adapt1_bandlimits1,'String'));
XStimParams.bandlimits(2) = str2num(get(H.Adapt1_bandlimits2,'String'));

% set numreps
XStimParams.numreps = str2num(get(H.Adapt1_numreps,'String'));

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
if(exist1('H.Adapt1fig') & get(H.Adapt1_inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.Adapt1_inc_testnum,'Value',0);
end
set(H.Adapt1_data_FN,'String', FN.data);

% update display of stim filename
stim_type = get(H.Adapt1_stim_type,'String');
stim_val = get(H.Adapt1_stim_type,'Value');
set(H.stim_type,'Value',stim_val);
XStimParams.stim_type = deblank(stim_type(stim_val,:));

clear stim_type

set(H.stim_filename,'Enable','off');
set(H.Adapt1_stimFN,'Enable','off');
set(H.Adapt1_freqtxt,'visible','off');
set(H.Adapt1_freq,'visible','off');
set(H.Adapt1_bandlimitstxt,'visible','off');
set(H.Adapt1_bandlimits1,'visible','off');
set(H.Adapt1_bandlimits2,'visible','off');

switch get(H.stim_type,'Value')
    case {1 2}
        set(H.Adapt1_freqtxt,'visible','on');
        set(H.Adapt1_freq,'visible','on');
        XStimParams.curr_freq = str2num(get(H.Adapt1_freq,'string'));
    case 6
        set(H.Adapt1_bandlimitstxt,'visible','on');
        set(H.Adapt1_bandlimits1,'visible','on');
        set(H.Adapt1_bandlimits2,'visible','on');
    case 8
    case 10
    otherwise
    set(H.stim_type,'Value',8);
    disp('Unsupported stimulus type for Adapt1 test; reset to BROADBAND');
end

% check out modulation parameters for stim1
XStimParams.mod_depth(1) = str2num(get(H.Adapt1_mod_depth(1),'String'));
mod_type = get(H.Adapt1_mod_type(1),'String');
mod_num = get(H.Adapt1_mod_type(1),'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Sq wave'
        set(H.Adapt1_mod_txt(1),'String', 'Freq:');
        set(H.Adapt1_mod_freq(1),'Visible','on');
        set(H.Adapt1_mod_phase(1),'Visible','on');
        set(H.Adapt1_mod_txtA(1),'Visible','on');
        set(H.Adapt1_mod_depth(1),'Visible','on');
        set(H.Adapt1_mod_txtB(1),'Visible','on');
        set(H.Adapt1_mod_pb(1),'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Adapt1_mod_freq(1),'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Adapt1_mod_phase(1),'String'));
    case 'Tone'
        set(H.Adapt1_mod_txt(1),'String', 'Freq:');
        set(H.Adapt1_mod_freq(1),'Visible','on');
        set(H.Adapt1_mod_phase(1),'Visible','on');
        set(H.Adapt1_mod_txtA(1),'Visible','on');
        set(H.Adapt1_mod_depth(1),'Visible','on');
        set(H.Adapt1_mod_txtB(1),'Visible','on');
        set(H.Adapt1_mod_pb(1),'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Adapt1_mod_freq(1),'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Adapt1_mod_phase(1),'String'));
    case 'File'
        set(H.Adapt1_mod_pb(1),'Visible','on');
        set(H.Adapt1_mod_freq(1),'Visible','off');
        set(H.Adapt1_mod_txtA(1),'Visible','on');
        set(H.Adapt1_mod_depth(1),'Visible','on');
        set(H.Adapt1_mod_phase(1),'Visible','off');
        set(H.Adapt1_mod_txtB(1),'Visible','off');
        if get(H.Adapt1_mod_pb(1),'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.Adapt1_mod_pb(1),'Value',0);
        end
        set(H.Adapt1_mod_txt(1),'String', FN.mod);
        
    case 'LP Noise'
        set(H.Adapt1_mod_txt(1),'String', 'CutOff Freq:');
        set(H.Adapt1_mod_freq(1),'Visible','on');
        set(H.Adapt1_mod_txtA(1),'Visible','on');
        set(H.Adapt1_mod_depth(1),'Visible','on');
        set(H.Adapt1_mod_phase(1),'Visible','off');
        set(H.Adapt1_mod_pb(1),'Visible','off');
        set(H.Adapt1_mod_txtB(1),'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Adapt1_mod_freq(1),'String'));
        
    case 'None'
        set(H.Adapt1_mod_txt(1),'String', 'none');
        set(H.Adapt1_mod_freq(1),'Visible','off');
        set(H.Adapt1_mod_phase(1),'Visible','off');
        set(H.Adapt1_mod_pb(1),'Visible','off');
        set(H.Adapt1_mod_txtB(1),'Visible','off');
        set(H.Adapt1_mod_txtA(1),'Visible','off');
        set(H.Adapt1_mod_depth(1),'Visible','off');
        
    otherwise
end

% check out modulation parameters for stim2
XStimParams.mod_depth(2) = str2num(get(H.Adapt1_mod_depth(2),'String'));
mod_type = get(H.Adapt1_mod_type(2),'String');
mod_num = get(H.Adapt1_mod_type(2),'Value');
XStimParams.mod_type2 = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type2
    case 'Sq wave'
        set(H.Adapt1_mod_txt(2),'String', 'Freq:');
        set(H.Adapt1_mod_freq(2),'Visible','on');
        set(H.Adapt1_mod_phase(2),'Visible','on');
        set(H.Adapt1_mod_txtA(2),'Visible','on');
        set(H.Adapt1_mod_depth(2),'Visible','on');
        set(H.Adapt1_mod_txtB(2),'Visible','on');
        set(H.Adapt1_mod_pb(2),'Visible','off');
        XStimParams.mod_freq(2) = str2num(get(H.Adapt1_mod_freq(2),'String'));
        XStimParams.mod_phase(2) = str2num(get(H.Adapt1_mod_phase(2),'String'));
    case 'Tone'
        set(H.Adapt1_mod_txt(2),'String', 'Freq:');
        set(H.Adapt1_mod_freq(2),'Visible','on');
        set(H.Adapt1_mod_phase(2),'Visible','on');
        set(H.Adapt1_mod_txtA(2),'Visible','on');
        set(H.Adapt1_mod_depth(2),'Visible','on');
        set(H.Adapt1_mod_txtB(2),'Visible','on');
        set(H.Adapt1_mod_pb(2),'Visible','off');
        XStimParams.mod_freq(2) = str2num(get(H.Adapt1_mod_freq(2),'String'));
        XStimParams.mod_phase(2) = str2num(get(H.Adapt1_mod_phase(2),'String'));
    case 'File'
        set(H.Adapt1_mod_pb(2),'Visible','on');
        set(H.Adapt1_mod_freq(2),'Visible','off');
        set(H.Adapt1_mod_txtA(2),'Visible','on');
        set(H.Adapt1_mod_depth(2),'Visible','on');
        set(H.Adapt1_mod_phase(2),'Visible','off');
        set(H.Adapt1_mod_txtB(2),'Visible','off');
        if get(H.Adapt1_mod_pb(2),'Value')
            [mod_from_file, FN.mod2, FN.mod_path2] = Get_Stimulus_File(XStimParams, FN.mod_path2);
            set(H.Adapt1_mod_pb(2),'Value',0);
        end
        set(H.Adapt1_mod_txt(2),'String', FN.mod2);
        
    case 'LP Noise'
        set(H.Adapt1_mod_txt(2),'String', 'CutOff Freq:');
        set(H.Adapt1_mod_freq(2),'Visible','on');
        set(H.Adapt1_mod_txtA(2),'Visible','on');
        set(H.Adapt1_mod_depth(2),'Visible','on');
        set(H.Adapt1_mod_phase(2),'Visible','off');
        set(H.Adapt1_mod_pb(2),'Visible','off');
        set(H.Adapt1_mod_txtB(2),'Visible','off');
        XStimParams.mod_freq(2) = str2num(get(H.Adapt1_mod_freq(2),'String'));
        
    case 'None'
        set(H.Adapt1_mod_txt(2),'String', 'none');
        set(H.Adapt1_mod_freq(2),'Visible','off');
        set(H.Adapt1_mod_phase(2),'Visible','off');
        set(H.Adapt1_mod_pb(2),'Visible','off');
        set(H.Adapt1_mod_txtB(2),'Visible','off');
        set(H.Adapt1_mod_txtA(2),'Visible','off');
        set(H.Adapt1_mod_depth(2),'Visible','off');
        
    otherwise
end

update_XStim;

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])