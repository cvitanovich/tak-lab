% SetInfo_space
% combines all the other 'set...' scripts into one

% uses .STD HRIRs, therefore make sure earphone filters are also selected
   XStimParams.ephone_flag = 1;
   set(H.ephoneuseit,'Value',1);

 
if get(H.space_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.space_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.space_DynamicRangepb,'BackgroundColor','red');
    set(H.space_DynamicRangeTxt,'visible','on');
else
    XStimParams.HiDynamicRange = 0;
    set(H.space_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.space_DynamicRangepb,'BackgroundColor','yellow');
    set(H.space_DynamicRangeTxt,'visible','off');
end
   
if get(H.space_buildOnly,'value') | exist1('H.composite_fig')
    XStimParams.buildOnly = 1;
else
    XStimParams.buildOnly = 0;
end
set(H.space_buildOnly,'value',XStimParams.buildOnly);
   
if get(H.spaceflag,'Value') & ~XStimParams.space_flag
    XStimParams.space_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
    XStimParams.ABLalone_flag = 0;
elseif get(H.space_ABLaloneflag,'Value') & ~XStimParams.ABLalone_flag
    XStimParams.space_flag = 0;
    XStimParams.ABLalone_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
elseif get(H.space_ildaloneflag,'Value') & ~XStimParams.ildalone_flag
    XStimParams.space_flag = 0;
    XStimParams.ildalone_flag = 1;
    XStimParams.itdalone_flag = 0;
    XStimParams.ABLalone_flag = 0;
elseif get(H.space_itdaloneflag,'Value') & ~XStimParams.itdalone_flag
    XStimParams.space_flag = 0;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 1;
    XStimParams.ABLalone_flag = 0;
end
set(H.spaceflag,'Value',XStimParams.space_flag);
set(H.space_ildaloneflag,'value',XStimParams.ildalone_flag);
set(H.space_itdaloneflag,'value',XStimParams.itdalone_flag);
set(H.space_ABLaloneflag,'value',XStimParams.ABLalone_flag);

% control visibility of file picking for ILA, ITA or fully-cued
if get(H.spaceflag,'Value')
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');   
    set(H.space_ildalonefilehdr,'Enable','off');
    set(H.space_ildalonefile,'Enable','off');
    set(H.space_ildalonepb,'Enable','off');  
    set(H.space_itdalonefilehdr,'Enable','off');
    set(H.space_itdalonefile,'Enable','off');
    set(H.space_itdalonepb,'Enable','off');  
    set(H.space_ITDtxt,'BackgroundColor',[.8 .8 .8]);
    set(H.space_ILDtxt,'BackgroundColor',[.8 .8 .8]);
    set(H.space_plotwts,'Enable','off');
elseif get(H.space_ABLaloneflag,'Value')
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');   
    set(H.space_ildalonefilehdr,'Enable','off');
    set(H.space_ildalonefile,'Enable','off');
    set(H.space_ildalonepb,'Enable','off');  
    set(H.space_itdalonefilehdr,'Enable','off');
    set(H.space_itdalonefile,'Enable','off');
    set(H.space_itdalonepb,'Enable','off');  
    set(H.space_ITDtxt,'BackgroundColor','green');
    set(H.space_ILDtxt,'BackgroundColor','green');
    set(H.space_plotwts,'Enable','off');
elseif get(H.space_ildaloneflag,'Value')
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');
    set(H.space_ildalonefilehdr,'Enable','on');
    set(H.space_ildalonefile,'Enable','on');
    set(H.space_ildalonepb,'Enable','on');
    set(H.space_itdalonefilehdr,'Enable','off');
    set(H.space_itdalonefile,'Enable','off');
    set(H.space_itdalonepb,'Enable','off');  
    set(H.space_ITDtxt,'BackgroundColor','green');
    set(H.space_ILDtxt,'BackgroundColor',[.8 .8 .8]);
    set(H.space_plotwts,'Enable','on');
elseif get(H.space_itdaloneflag,'Value')
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');
    set(H.space_ildalonefilehdr,'Enable','off');
    set(H.space_ildalonefile,'Enable','off');
    set(H.space_ildalonepb,'Enable','off');
    set(H.space_itdalonefilehdr,'Enable','on');
    set(H.space_itdalonefile,'Enable','on');
    set(H.space_itdalonepb,'Enable','on');  
    set(H.space_ITDtxt,'BackgroundColor',[.8 .8 .8]);
    set(H.space_ILDtxt,'BackgroundColor','green');
    set(H.space_plotwts,'Enable','off');
end

% check for ILA pushbutton (file-picker)
if get(H.space_ildalonepb,'Value')
    [FN.ildalone_std,FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select ILD-alone HRTF File');
    if(FN.ILA_path ~= 0)
        set(H.space_ildalonefile,'String',FN.ildalone_std);
    end
    set(H.space_ildalonepb,'Value',0);
    FN.HRTFfiletype(2,2) = testHRTFfiletype(FN.ILA_path, FN.ildalone_std);
end

% check for ITA pushbutton (file-picker)
if get(H.space_itdalonepb,'Value')
    [FN.itdalone_std,FN.ITA_path] = uigetfile([FN.ITA_path '*.*'],'Select ITD-alone HRTF File');
    if(FN.ITA_path ~= 0)
        set(H.space_itdalonefile,'String',FN.itdalone_std);
    end
    set(H.space_itdalonepb,'Value',0);
    FN.HRTFfiletype(3,2) = testHRTFfiletype(FN.ITA_path, FN.itdalone_std);
end

% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
    [FN.space_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File (*.std preferred)');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space_std);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
end

% set reference locations for ILA
if get(H.spaceflag,'Value') | get(H.space_ABLaloneflag,'Value') | get(H.space_itdaloneflag,'Value')
    XStimParams.az = str2num(get(H.space_az,'String'));
    XStimParams.el = str2num(get(H.space_el,'String'));
end

% reset curr_ILD and curr_ITD
XStimParams.curr_ITD = str2num(get(H.space_ITD,'string'));
XStimParams.curr_ABL = round(str2num(get(H.space_ABL,'string')));
set(H.space_ABL,'string',num2str(XStimParams.curr_ABL));

% bandwidth & limits
if get(H.stim_type,'Value') ~= 9			% not from file
    XStimParams.bandwidth = str2num(get(H.space_bandwidth,'String'));
    XStimParams.bandlimits(1) = str2num(get(H.space_bandlimits1,'String'));
    XStimParams.bandlimits(2) = str2num(get(H.space_bandlimits2,'String'));
end

% set numreps
XStimParams.numreps = str2num(get(H.space_numreps,'String'));

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
if(exist1('H.spacefig') & get(H.space_inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.space_inc_testnum,'Value',0);
end
set(H.spacedata_FN,'String', FN.data);

% update display of stim filename
stim_type = get(H.space_stim_type,'String');
stim_val = get(H.space_stim_type,'Value');
set(H.stim_type,'Value',stim_val);
XStimParams.stim_type = deblank(stim_type(stim_val,:));
if get(H.space_stimFNpb,'value')
   [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
end

clear stim_type
if get(H.stim_type,'Value') ~= 7
    set(H.space_NBstackpb,'Visible','off');
    set(H.space_bandtxt(1),'Visible','off');
    set(H.space_bandwidth,'Visible','off');
    set(H.space_bandtxt(2),'Visible','off');
    set(H.space_bandlimits1,'Visible','off');
    set(H.space_bandlimits2,'Visible','off');
    set(H.space_bandtxt(3),'Visible','on');
    set(H.space_curr_freq,'Visible','on');
else
    set(H.space_NBstackpb,'Visible','on');
    set(H.space_bandtxt(1),'Visible','on');
    set(H.space_bandwidth,'Visible','on');
    set(H.space_bandtxt(2),'Visible','on');
    set(H.space_bandlimits1,'Visible','on');
    set(H.space_bandlimits2,'Visible','on');
    set(H.space_bandtxt(3),'Visible','off');
    set(H.space_curr_freq,'Visible','off');
end
XStimParams.curr_freq = str2num(get(H.space_curr_freq,'string'));

if get(H.stim_type,'Value') == 9								% from file
    set(H.space_stimFN,'String',FN.stim);
    set(H.stim_filename,'String',FN.stim);
    set(H.stim_filename,'Enable','on');
    set(H.space_stimFN,'Enable','on');
    D = dir([FN.stim_path FN.stim]);
    XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.space_DUR, 'String', num2str(XStimParams.curr_stimdur));
else
    set(H.stim_filename,'Enable','off');
    set(H.space_stimFN,'Enable','off');
end

% update stimulus duration
XStimParams.curr_stimdur = str2num(get(H.space_DUR, 'String'));

% check out modulation parameters for stim
XStimParams.mod_depth(1) = str2num(get(H.spacetestmod_depth,'String'));
mod_type = get(H.spacetestmod_type,'String');
mod_num = get(H.spacetestmod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Sq wave'
        set(H.spacetestmod_txt,'String', 'Freq:');
        set(H.spacetestmod_freq,'Visible','on');
        set(H.spacetestmod_phase,'Visible','on');
        set(H.spacetestmod_txtA,'Visible','on');
        set(H.spacetestmod_depth,'Visible','on');
        set(H.spacetestmod_txtB,'Visible','on');
        set(H.spacetestmod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.spacetestmod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.spacetestmod_phase,'String'));
    case 'Tone'
        set(H.spacetestmod_txt,'String', 'Freq:');
        set(H.spacetestmod_freq,'Visible','on');
        set(H.spacetestmod_phase,'Visible','on');
        set(H.spacetestmod_txtA,'Visible','on');
        set(H.spacetestmod_depth,'Visible','on');
        set(H.spacetestmod_txtB,'Visible','on');
        set(H.spacetestmod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.spacetestmod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.spacetestmod_phase,'String'));
    case 'File'
        set(H.spacetestmod_pb,'Visible','on');
        set(H.spacetestmod_freq,'Visible','off');
        set(H.spacetestmod_txtA,'Visible','on');
        set(H.spacetestmod_depth,'Visible','on');
        set(H.spacetestmod_phase,'Visible','off');
        set(H.spacetestmod_txtB,'Visible','off');
        if get(H.spacetestmod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.spacetestmod_pb,'Value',0);
        end
        set(H.spacetestmod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.spacetestmod_txt,'String', 'CutOff Freq:');
        set(H.spacetestmod_freq,'Visible','on');
        set(H.spacetestmod_txtA,'Visible','on');
        set(H.spacetestmod_depth,'Visible','on');
        set(H.spacetestmod_phase,'Visible','off');
        set(H.spacetestmod_pb,'Visible','off');
        set(H.spacetestmod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.spacetestmod_freq,'String'));
        
    case 'None'
        set(H.spacetestmod_txt,'String', 'no mod  ');
        set(H.spacetestmod_freq,'Visible','off');
        set(H.spacetestmod_phase,'Visible','off');
        set(H.spacetestmod_pb,'Visible','off');
        set(H.spacetestmod_txtB,'Visible','off');
        set(H.spacetestmod_txtA,'Visible','off');
        set(H.spacetestmod_depth,'Visible','off');
        
    otherwise
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim