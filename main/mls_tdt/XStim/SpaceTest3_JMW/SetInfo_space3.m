% SetInfo_space3
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
    [FN.ildalone_eq,FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select ILD-alone HRTF File');
    if(FN.ILA_path ~= 0)
        set(H.ildalonefile,'String',FN.ildalone_eq);
    end
    set(H.ildalonepb,'Value',0);
    FN.HRTFfiletype(2,1) = testHRTFfiletype(FN.ILA_path, FN.ildalone_eq);
end

% check for ITA pushbutton (file-picker)
if get(H.itdalonepb,'Value')
    [FN.itdalone_eq,FN.ITA_path] = uigetfile([FN.ITA_path '*.*'],'Select ITD-alone HRTF File');
    if(FN.ITA_path ~= 0)
        set(H.itdalonefile,'String',FN.itdalone_eq);
    end
    set(H.itdalonepb,'Value',0);
    FN.HRTFfiletype(3,1) = testHRTFfiletype(FN.ITA_path, FN.itdalone_eq);
end

% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
    [FN.space_eq,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select Fully-cued HRTF File');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space_eq);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);
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
if(exist1('H.spacefig3') & get(H.inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.inc_testnum,'Value',0);
end

if get(H.space3_stimFNpb,'Value')
    [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
    set(H.space3_stimFNpb,'value',0);
end

% update display of stim filename
if get(H.space3_stim_type,'Value') == 9								% from file
    set(H.space3_stimFN,'String',FN.stim);
    set(H.stim_filename,'String',FN.stim);
    set(H.stim_filename,'Enable','on');
    set(H.space3_stimFN,'Enable','on');
    D = dir([FN.stim_path FN.stim]);
    XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.DUR, 'String', num2str(XStimParams.curr_stimdur));
else
    set(H.stim_filename,'Enable','off');
    set(H.space3_stimFN,'Enable','off');
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