% SetInfo_PlayFile
% combines all the other 'set...' scripts into one

% uses .STD HRIRs, therefore make sure earphone filters are also selected
   XStimParams.ephone_flag = 1;
   set(H.ephoneuseit,'Value',1);

 
if get(H.PlayFile_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.PlayFile_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.PlayFile_DynamicRangepb,'BackgroundColor','red');
    set(H.PlayFile_DynamicRangeTxt,'visible','on');
else
    XStimParams.HiDynamicRange = 0;
    set(H.PlayFile_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.PlayFile_DynamicRangepb,'BackgroundColor','yellow');
    set(H.PlayFile_DynamicRangeTxt,'visible','off');
end
 
if get(H.PlayFile_NoScalingpb,'value')
    set(H.PlayFile_NoScalingpb,'string','push for Normal Stim Scaling');
    set(H.PlayFile_NoScalingpb,'BackgroundColor','red');
else
    set(H.PlayFile_NoScalingpb,'string','push for No Stim Scaling');
    set(H.PlayFile_NoScalingpb,'BackgroundColor','yellow');
end




% set reference locations for ILA

% reset curr_ILD and curr_ITD
XStimParams.curr_ITD = str2num(get(H.PlayFile_ITD,'string'));
XStimParams.curr_ABL = str2num(get(H.PlayFile_ABL,'string'));
XStimParams.curr_ILD = str2num(get(H.PlayFile_ILD,'string'));

% bandwidth & limits
if get(H.stim_type,'Value') ~= 9			% not from file
    XStimParams.bandwidth = str2num(get(H.PlayFile_bandwidth,'String'));
    XStimParams.bandlimits(1) = str2num(get(H.PlayFile_bandlimits1,'String'));
    XStimParams.bandlimits(2) = str2num(get(H.PlayFile_bandlimits2,'String'));
end

% set numreps
XStimParams.numreps = str2num(get(H.PlayFile_numreps,'String'));

%Set PlayFile PickerFlag
XStimParams.picker_flag = 0;

% increment test number
if(exist1('H.PlayFilefig') & get(H.PlayFile_inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.PlayFile_inc_testnum,'Value',0);
end
set(H.PlayFiledata_FN,'String', FN.data);

% update display of stim filename
stim_type = get(H.PlayFile_stim_type,'String');
stim_val = get(H.PlayFile_stim_type,'Value');
set(H.stim_type,'Value',stim_val);
XStimParams.stim_type = deblank(stim_type(stim_val,:));
if get(H.PlayFile_stimFNpb,'value')
   [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
end

clear stim_type

if get(H.stim_type,'Value') == 9								% from file
    set(H.PlayFile_stimFN,'String',FN.stim);
    set(H.stim_filename,'String',FN.stim);
    set(H.stim_filename,'Enable','on');
    set(H.PlayFile_stimFN,'Enable','on');
    D = dir([FN.stim_path FN.stim]);
    XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.PlayFile_DUR, 'String', num2str(XStimParams.curr_stimdur));
else
    set(H.stim_filename,'Enable','off');
    set(H.PlayFile_stimFN,'Enable','off');
end

% update stimulus duration
XStimParams.curr_stimdur = str2num(get(H.PlayFile_DUR, 'String'));

% check out modulation parameters for stim
XStimParams.mod_depth(1) = str2num(get(H.PlayFiletestmod_depth,'String'));
mod_type = get(H.PlayFiletestmod_type,'String');
mod_num = get(H.PlayFiletestmod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Sq wave'
        set(H.PlayFiletestmod_txt,'String', 'Freq:');
        set(H.PlayFiletestmod_freq,'Visible','on');
        set(H.PlayFiletestmod_phase,'Visible','on');
        set(H.PlayFiletestmod_txtA,'Visible','on');
        set(H.PlayFiletestmod_depth,'Visible','on');
        set(H.PlayFiletestmod_txtB,'Visible','on');
        set(H.PlayFiletestmod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.PlayFiletestmod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.PlayFiletestmod_phase,'String'));
    case 'Tone'
        set(H.PlayFiletestmod_txt,'String', 'Freq:');
        set(H.PlayFiletestmod_freq,'Visible','on');
        set(H.PlayFiletestmod_phase,'Visible','on');
        set(H.PlayFiletestmod_txtA,'Visible','on');
        set(H.PlayFiletestmod_depth,'Visible','on');
        set(H.PlayFiletestmod_txtB,'Visible','on');
        set(H.PlayFiletestmod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.PlayFiletestmod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.PlayFiletestmod_phase,'String'));
    case 'File'
        set(H.PlayFiletestmod_pb,'Visible','on');
        set(H.PlayFiletestmod_freq,'Visible','off');
        set(H.PlayFiletestmod_txtA,'Visible','on');
        set(H.PlayFiletestmod_depth,'Visible','on');
        set(H.PlayFiletestmod_phase,'Visible','off');
        set(H.PlayFiletestmod_txtB,'Visible','off');
        if get(H.PlayFiletestmod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.PlayFiletestmod_pb,'Value',0);
        end
        set(H.PlayFiletestmod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.PlayFiletestmod_txt,'String', 'CutOff Freq:');
        set(H.PlayFiletestmod_freq,'Visible','on');
        set(H.PlayFiletestmod_txtA,'Visible','on');
        set(H.PlayFiletestmod_depth,'Visible','on');
        set(H.PlayFiletestmod_phase,'Visible','off');
        set(H.PlayFiletestmod_pb,'Visible','off');
        set(H.PlayFiletestmod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.PlayFiletestmod_freq,'String'));
        
    case 'None'
        set(H.PlayFiletestmod_txt,'String', 'no mod  ');
        set(H.PlayFiletestmod_freq,'Visible','off');
        set(H.PlayFiletestmod_phase,'Visible','off');
        set(H.PlayFiletestmod_pb,'Visible','off');
        set(H.PlayFiletestmod_txtB,'Visible','off');
        set(H.PlayFiletestmod_txtA,'Visible','off');
        set(H.PlayFiletestmod_depth,'Visible','off');
        
    otherwise
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])