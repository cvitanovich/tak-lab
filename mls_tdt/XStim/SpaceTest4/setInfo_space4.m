% SetInfo_space4
% combines all the other 'set...' scripts into one

% uses .eq HRIRs, therefore make sure earphone filters are NOT also selected
XStimParams.ephone_flag = 0;
set(H.ephoneuseit,'Value',0);

% space4_LimitXCpb
% not functional yet...
if 0 
    if get(H.space4_LimitXCpb,'value')
        XStimParams.test_type = 'Space4_limitXC';
        set(H.space4_LimitXCpb,'string','push for Normal test');
        set(H.space4_LimitXCpb,'BackgroundColor','red');
        set(H.space4_LimitXCTxt,'string','XC b/wn freqs LIMITED');
    else
        XStimParams.test_type = 'Space4';
        set(H.space4_LimitXCpb,'string','push for XClimit test');
        set(H.space4_LimitXCpb,'BackgroundColor','yellow');
        set(H.space4_LimitXCTxt,'string','XC b/wn freqs NOT limited');
    end
else
    set(H.space4_LimitXCpb,'visible','off')
    set(H.space4_LimitXCTxt,'visible','off')
    XStimParams.test_type = 'Space4';
end

if get(H.space4_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.space4_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.space4_DynamicRangepb,'BackgroundColor','red');
    set(H.space4_DynamicRangeTxt,'visible','on');
else
    XStimParams.HiDynamicRange = 0;
    set(H.space4_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.space4_DynamicRangepb,'BackgroundColor','yellow');
    set(H.space4_DynamicRangeTxt,'visible','off');
end


XStimParams.buildOnly = 0;

XStimParams.space_flag = 1;
XStimParams.ildalone_flag = 0;
XStimParams.itdalone_flag = 0;
XStimParams.ABLalone_flag = 0;

set(H.space4flag,'Value',XStimParams.space_flag);

XStimParams.curr_ABL = str2num(get(H.space4_ABL,'string'));

% check for fully-cued pushbutton (file-picker)
if get(H.space4filepb,'Value')
    [FN.space_eq,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File (*.eq preferred)');
    if(FN.space_path ~= 0)
        set(H.space4file,'String',FN.space_eq);
    end
    set(H.space4filepb,'Value',0);
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_eq);
end

% set numreps
XStimParams.numreps = str2num(get(H.space4_numreps,'String'));

%Set space4 PickerFlag
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
if(exist1('H.space4fig') & get(H.space4_inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.space4_inc_testnum,'Value',0);
end
set(H.space4data_FN,'String', FN.data);

% update display of stim filename
stim_type = get(H.space4_stim_type,'String');
stim_val = get(H.space4_stim_type,'Value');
set(H.stim_type,'Value',stim_val);
XStimParams.stim_type = deblank(stim_type(stim_val,:));
if get(H.space4_stimFNpb,'value')
   [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
end

%    'Tone|Gammatone||||NarrowBand Noise||BroadBand Noise|File|BBN - LNN||',...
clear stim_type
if get(H.stim_type,'Value') == 9								% from file
    set(H.space4_stimFN,'String',FN.stim);
    set(H.stim_filename,'String',FN.stim);
    set(H.stim_filename,'Enable','on');
    set(H.space4_stimFN,'Enable','on');
    D = dir([FN.stim_path FN.stim]);
    XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.space4_DUR, 'String', num2str(XStimParams.curr_stimdur));
    set(H.space4_FNseriestxt,'Visible','on')
else
    set(H.stim_filename,'Enable','off');
    set(H.space4_stimFN,'Enable','off');
    set(H.space4_FNseriestxt,'Visible','off')
end


XStimParams.curr_freq = str2num(get(H.space4_curr_freq,'string'));

if get(H.stim_type,'Value')<=6
    set(H.space4_curr_freqtxt, 'visible','on');
    set(H.space4_curr_freq, 'visible','on');
else
    set(H.space4_curr_freqtxt, 'visible','off');
    set(H.space4_curr_freq, 'visible','off');    
end

% check out modulation parameters for stim
XStimParams.mod_depth(1) = str2num(get(H.space4testmod_depth,'String'));
mod_type = get(H.space4testmod_type,'String');
mod_num = get(H.space4testmod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Sq wave'
        set(H.space4testmod_txt,'String', 'Freq:');
        set(H.space4testmod_freq,'Visible','on');
        set(H.space4testmod_phase,'Visible','on');
        set(H.space4testmod_txtA,'Visible','on');
        set(H.space4testmod_depth,'Visible','on');
        set(H.space4testmod_txtB,'Visible','on');
        set(H.space4testmod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.space4testmod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.space4testmod_phase,'String'));
    case 'Tone'
        set(H.space4testmod_txt,'String', 'Freq:');
        set(H.space4testmod_freq,'Visible','on');
        set(H.space4testmod_phase,'Visible','on');
        set(H.space4testmod_txtA,'Visible','on');
        set(H.space4testmod_depth,'Visible','on');
        set(H.space4testmod_txtB,'Visible','on');
        set(H.space4testmod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.space4testmod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.space4testmod_phase,'String'));
    case 'File'
        set(H.space4testmod_pb,'Visible','on');
        set(H.space4testmod_freq,'Visible','off');
        set(H.space4testmod_txtA,'Visible','on');
        set(H.space4testmod_depth,'Visible','on');
        set(H.space4testmod_phase,'Visible','off');
        set(H.space4testmod_txtB,'Visible','off');
        if get(H.space4testmod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.space4testmod_pb,'Value',0);
        end
        set(H.space4testmod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.space4testmod_txt,'String', 'CutOff Freq:');
        set(H.space4testmod_freq,'Visible','on');
        set(H.space4testmod_txtA,'Visible','on');
        set(H.space4testmod_depth,'Visible','on');
        set(H.space4testmod_phase,'Visible','off');
        set(H.space4testmod_pb,'Visible','off');
        set(H.space4testmod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.space4testmod_freq,'String'));
        
    case 'None'
        set(H.space4testmod_txt,'String', 'no mod  ');
        set(H.space4testmod_freq,'Visible','off');
        set(H.space4testmod_phase,'Visible','off');
        set(H.space4testmod_pb,'Visible','off');
        set(H.space4testmod_txtB,'Visible','off');
        set(H.space4testmod_txtA,'Visible','off');
        set(H.space4testmod_depth,'Visible','off');
        
    otherwise
end

XStimParams.mu = str2num(get(H.space4_mu,'string'));
XStimParams.sigma = str2num(get(H.space4_sigma,'string'));
XStimParams.epoch_duration(1) = str2num(get(H.space4_epoch_duration(1),'string'));
XStimParams.focalProb = str2num(get(H.space4_focalProb,'string'));

% update stimulus duration
if get(H.stim_type,'Value') ~= 9
    XStimParams.curr_stimdur = str2num(get(H.space4_DUR, 'String'));
    nEpochs = round(XStimParams.curr_stimdur / XStimParams.epoch_duration(1));
    XStimParams.curr_stimdur = nEpochs * XStimParams.epoch_duration(1);
    if nEpochs ~= length(XStimParams.mu)
        set(H.space4_mu,'backgroundcolor','red')
        Htemp = errordlg('change mu or stimdur','nEpochs does not equal length(mu)');
        set(Htemp,'position',[400 500 200 100]);
    else
        set(H.space4_mu,'backgroundcolor',[ 0.6275    0.7216    0.7843])
    end
    if nEpochs ~= length(XStimParams.sigma)
        set(H.space4_sigma,'backgroundcolor','red')
        Htemp = errordlg('Change sigma or stimdur','nEpochs does not equal length(sigma)');
        set(Htemp,'position',[400 500 200 100]);
    else
        set(H.space4_sigma,'backgroundcolor',[ 0.6275    0.7216    0.7843])
    end
else
    nEpochs = 1;
end

set(H.space4_DUR,'string',num2str(XStimParams.curr_stimdur));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim