% SetInfo_Two_source
% combines all the other 'set...' scripts into one

%check ILDAloneFlag
% XStimParams.space_flag == 1 for ILDalone
%						 == 0 for fully cued

if get(H.Twosrc_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.Twosrc_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.Twosrc_DynamicRangepb,'BackgroundColor','red');
    set(H.Twosrc_DynamicRangeTxt,'visible','on');
else
    XStimParams.HiDynamicRange = 0;
    set(H.Twosrc_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.Twosrc_DynamicRangepb,'BackgroundColor','yellow');
    set(H.Twosrc_DynamicRangeTxt,'visible','off');
end

if get(H.spaceflag,'Value') & (XStimParams.space_flag == 0)
    XStimParams.space_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
    XStimParams.ABLalone_flag = 0;
elseif get(H.ABLaloneflag,'Value') & (XStimParams.ABLalone_flag == 0)
    XStimParams.space_flag = 0;
    XStimParams.ABLalone_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
elseif get(H.ildaloneflag,'Value') & (XStimParams.ildalone_flag == 0)
    XStimParams.space_flag = 0;
    XStimParams.ildalone_flag = 1;
    XStimParams.itdalone_flag = 0;
    XStimParams.ABLalone_flag = 0;
elseif get(H.itdaloneflag,'Value') & (XStimParams.itdalone_flag == 0)
    XStimParams.space_flag = 0;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 1;
    XStimParams.ABLalone_flag = 0;
end
set(H.spaceflag,'Value',XStimParams.space_flag);
set(H.ildaloneflag,'value',XStimParams.ildalone_flag);
set(H.itdaloneflag,'value',XStimParams.itdalone_flag);
set(H.ABLaloneflag,'value',XStimParams.ABLalone_flag);

% control visibility of file picking for ILA, ITA or fully-cued
% control visibility of file picking for ILA, ITA or fully-cued
if XStimParams.space_flag
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');   
    set(H.ildalonefilehdr,'Enable','off');
    set(H.ildalonefile,'Enable','off');
    set(H.ildalonepb,'Enable','off');  
    set(H.itdalonefilehdr,'Enable','off');
    set(H.itdalonefile,'Enable','off');
    set(H.itdalonepb,'Enable','off');  
    set(H.ITDtxt,'BackgroundColor',[.8 .8 .8]);
    set(H.ILDtxt,'BackgroundColor',[.8 .8 .8]);
elseif XStimParams.ABLalone_flag
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');   
    set(H.ildalonefilehdr,'Enable','off');
    set(H.ildalonefile,'Enable','off');
    set(H.ildalonepb,'Enable','off');  
    set(H.itdalonefilehdr,'Enable','off');
    set(H.itdalonefile,'Enable','off');
    set(H.itdalonepb,'Enable','off');  
    set(H.ITDtxt,'BackgroundColor','green');
    set(H.ILDtxt,'BackgroundColor','green');
elseif XStimParams.ildalone_flag
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');
    set(H.ildalonefilehdr,'Enable','on');
    set(H.ildalonefile,'Enable','on');
    set(H.ildalonepb,'Enable','on');
    set(H.itdalonefilehdr,'Enable','off');
    set(H.itdalonefile,'Enable','off');
    set(H.itdalonepb,'Enable','off');  
    set(H.ITDtxt,'BackgroundColor','green');
    set(H.ILDtxt,'BackgroundColor',[.8 .8 .8]);
elseif XStimParams.itdalone_flag
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');
    set(H.ildalonefilehdr,'Enable','off');
    set(H.ildalonefile,'Enable','off');
    set(H.ildalonepb,'Enable','off');
    set(H.itdalonefilehdr,'Enable','on');
    set(H.itdalonefile,'Enable','on');
    set(H.itdalonepb,'Enable','on');  
    set(H.ITDtxt,'BackgroundColor',[.8 .8 .8]);
    set(H.ILDtxt,'BackgroundColor','green');
end

% check for IA pushbutton (file-picker)
if get(H.ildalonepb,'Value')
    [FN.ildalone_std,FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select ILD-alone HRTF File');
    if(FN.ILA_path ~= 0)
        set(H.ildalonefile,'String',FN.ildalone_std);
    end
    set(H.ildalonepb,'Value',0);
    FN.HRTFfiletype(2,2) = testHRTFfiletype(FN.ILA_path, FN.ildalone_std);
end

% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
    [FN.space_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space_std);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
end
% check for ITA pushbutton (file-picker)
if get(H.itdalonepb,'Value')
    [FN.itdalone,FN.ITA_path] = uigetfile([FN.ITA_path '*.*'],'Select ITD-alone HRTF File');
    if(FN.ITA_path ~= 0)
        set(H.itdalonefile,'String',FN.itdalone_std);
    end
    set(H.itdalonepb,'Value',0);
    FN.HRTFfiletype(3,2) = testHRTFfiletype(FN.ITA_path, FN.itdalone_std);
end

% set reference locations
if get(H.spaceflag,'Value') | get(H.ABLaloneflag,'Value') | get(H.itdaloneflag,'Value')
    XStimParams.az = str2num(get(H.az,'String'));
    XStimParams.el = str2num(get(H.el,'String'));
end

% set numreps
XStimParams.numreps = str2num(get(H.numreps,'String'));

% set offsets
XStimParams.offset_el(1) = 0;
XStimParams.offset_az(1) = 0;
XStimParams.offset_el(2) = str2num(get(H.offset_el,'String'));
XStimParams.offset_az(2) = str2num(get(H.offset_az,'String'));
XStimParams.offset_el(3) = str2num(get(H.offset_el2,'String'));
XStimParams.offset_az(3) = str2num(get(H.offset_az2,'String'));

% set multiplication factor for each source
XStimParams.factor(1) = str2num(get(H.Twosrc_factor,'String'));
XStimParams.factor(2) = str2num(get(H.Twosrc_factor2,'String'));
XStimParams.factor(3) = str2num(get(H.Twosrc_factor3,'String'));

% set Use of third stim
XStimParams.ThreeStims = get(H.Use3Stims,'Value');
if XStimParams.ThreeStims
    set(H.Use3Stims_txt, 'BackgroundColor','yellow');
    set(H.Use3Stims, 'BackgroundColor','yellow');
    set(H.Stims_txt1, 'BackgroundColor','yellow');
    set(H.offset_el2, 'BackgroundColor','yellow');
    set(H.Stims_txt2, 'BackgroundColor','yellow');
    set(H.offset_az2, 'BackgroundColor','yellow');
else
    set(H.Use3Stims_txt, 'BackgroundColor',[.4 .7 .7]);
    set(H.Use3Stims, 'BackgroundColor',[.4 .7 .7]);
    set(H.Stims_txt1, 'BackgroundColor',[.4 .7 .7]);
    set(H.offset_el2, 'BackgroundColor',[.4 .7 .7]);
    set(H.Stims_txt2, 'BackgroundColor',[.4 .7 .7]);
    set(H.offset_az2, 'BackgroundColor',[.4 .7 .7]);
end    


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
if(exist1('H.Two_sourcefig') & get(H.inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.inc_testnum,'Value',0);
end

% update stim_type1
stim_type = get(H.Twosrc_stim_type,'String');
stim_val = get(H.Twosrc_stim_type,'Value');

if ~strcmp(XStimParams.stim_type, deblank(stim_type(stim_val,:))) | ...
        get(H.Twosrc_stim_filenamepb,'Value')		% stim_type has changed or pb pushed
    XStimParams.stim_type = deblank(stim_type(stim_val,:));
    set(H.stim_type,'Value',stim_val);
    if(strcmp('File',deblank(XStimParams.stim_type)))
        [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
        set(H.Twosrc_stim_filename,'String', FN.stim);
        set(H.stim_filename,'String', FN.stim);
    end
    %setinfo_Xstim;
    set(H.Twosrc_stim_filenamepb,'Value',0);
end

%control visibility of FN vs freq
set(H.Twosrc_stim_filenamepb,'Visible','on');
set(H.Twosrc_stim_filename,'Visible','on');
set(H.Twosrc_freq_txt,'Visible','off');
set(H.Twosrc_freq,'Visible','off');    

switch deblank(XStimParams.stim_type)
    case 'File'
        set(H.Twosrc_stim_filename,'Enable','on');
        set(H.stim_filename,'Enable','on');
    case {'Gammatone', 'Tone'}
        set(H.Twosrc_freq_txt,'Visible','on');
        set(H.Twosrc_freq,'Visible','on');    
        set(H.Twosrc_stim_filenamepb,'Visible','off');
        set(H.Twosrc_stim_filename,'Visible','off');
    otherwise
        set(H.Twosrc_stim_filename,'Enable','off');
        set(H.stim_filename,'Enable','off');
end

% update stim_type2
stim_type = get(H.Twosrc_stim_type2,'String');
stim_val = get(H.Twosrc_stim_type2,'Value');

if ~strcmp(XStimParams.stim_type2,deblank(stim_type(stim_val,:))) | ...
        get(H.Twosrc_stim_filename2pb,'Value')		% stim_type2 has changed
    XStimParams.stim_type2 = deblank(stim_type(stim_val,:));
    set(H.stim_type2,'Value',stim_val);
    if(strcmp('File',deblank(XStimParams.stim_type2)))
        [stim_from_file2, FN.stim2, FN.stim_path2] = Get_Stimulus_File(XStimParams, FN.stim_path2);
        set(H.Twosrc_stim_filename2,'String', FN.stim2);
        set(H.stim_filename2,'String', FN.stim2);
    end
    set(H.Twosrc_stim_filename2pb,'Value',0);
end
%control visibility of FN vs freq
set(H.Twosrc_stim_filename2pb,'Visible','on');
set(H.Twosrc_stim_filename2,'Visible','on');
set(H.Twosrc_freq2_txt,'Visible','off');
set(H.Twosrc_freq2,'Visible','off');    

switch deblank(XStimParams.stim_type2)
    case 'File'
        set(H.Twosrc_stim_filename2,'Enable','on');
        set(H.stim_filename2,'Enable','on');
    case {'Gammatone', 'Tone'}
        set(H.Twosrc_freq2_txt,'Visible','on');
        set(H.Twosrc_freq2,'Visible','on');    
        set(H.Twosrc_stim_filename2pb,'Visible','off');
        set(H.Twosrc_stim_filename2,'Visible','off');
    otherwise
        set(H.Twosrc_stim_filename2,'Enable','off');
        set(H.stim_filename2,'Enable','off');
end

% update stim_type3
stim_type = get(H.Twosrc_stim_type3,'String');
stim_val = get(H.Twosrc_stim_type3,'Value');

if ~strcmp(XStimParams.stim_type3, deblank(stim_type(stim_val,:))) | ...
        get(H.Twosrc_stim_filenamepb3,'Value')		% stim_type has changed or pb pushed
    XStimParams.stim_type3 = deblank(stim_type(stim_val,:));
    set(H.stim_type3,'Value',stim_val);
    if(strcmp('File',deblank(XStimParams.stim_type3)))
        [stim_from_file, FN.stim3, FN.stim_path3] = Get_Stimulus_File(XStimParams, FN.stim_path3);
        set(H.Twosrc_stim_filename3,'String', FN.stim3);
        set(H.stim_filename3,'String', FN.stim3);
    end
    set(H.Twosrc_stim_filenamepb3,'Value',0);
end
%control visibility of FN vs freq
set(H.Twosrc_stim_filenamepb3,'Visible','on');
set(H.Twosrc_stim_filename3,'Visible','on');
set(H.Twosrc_freq_txt3,'Visible','off');
set(H.Twosrc_freq3,'Visible','off');    

switch deblank(XStimParams.stim_type3)
    case 'File'
        set(H.Twosrc_stim_filename3,'Enable','on');
        set(H.stim_filename3,'Enable','on');
    case {'Gammatone', 'Tone'}
        set(H.Twosrc_freq_txt3,'Visible','on');
        set(H.Twosrc_freq3,'Visible','on');    
        set(H.Twosrc_stim_filenamepb3,'Visible','off');
        set(H.Twosrc_stim_filename3,'Visible','off');
    otherwise
        set(H.Twosrc_stim_filename3,'Enable','off');
        set(H.stim_filename3,'Enable','off');
end



% set dur based on larger of stim1, stim2 or stim3
set(H.DUR,'Enable','on');
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
temp = zeros(4,1);
temp(1) = XStimParams.curr_stimdur;
if (get(H.Twosrc_stim_type,'Value') == 9)       % from file
    set(H.DUR,'Enable','off');
    clear dir
    D = dir([FN.stim_path FN.stim]);
    temp(2) = round(1000 * (D.bytes/4) / TDT.Fs);
end
if (get(H.Twosrc_stim_type2,'Value') == 9)       % from file
    set(H.DUR,'Enable','off');
    clear dir
    D = dir([FN.stim_path2 FN.stim2]);
    temp(3) = round(1000 * (D.bytes/4) / TDT.Fs);
end
if (get(H.Twosrc_stim_type3,'Value') == 9) & XStimParams.ThreeStims      % from file
    set(H.DUR,'Enable','off');
    clear dir
    D = dir([FN.stim_path3 FN.stim3]);
    temp(4) = round(1000 * (D.bytes/4) / TDT.Fs);
end
XStimParams.curr_stimdur = max1(temp);
set(H.DUR,'string',num2str(XStimParams.curr_stimdur));


% bandwidths and bandlimits
% 1st stim
if ~strcmp(XStimParams.stim_type, 'File') 	% if not from file
    % if changed, re-set bandwidths
    if XStimParams.bandwidth(1) ~= str2num(get(H.bandwidth1,'String'))
        XStimParams.bandwidth(1) = str2num(get(H.bandwidth1,'String'));
    end
    % if changed, re-set bandlimits
    if (XStimParams.bandlimits(1) ~= str2num(get(H.bandlimits1,'String')) | ...
            XStimParams.bandlimits(2) ~= str2num(get(H.bandlimits2,'String')));  
        XStimParams.bandlimits(1) = str2num(get(H.bandlimits1,'String'));
        XStimParams.bandlimits(2) = str2num(get(H.bandlimits2,'String'));
    end
    set(H.bandwidth1,'Enable','on');
    set(H.bandlimits1,'Enable','on');
    set(H.bandlimits2,'Enable','on');
else
    set(H.bandwidth1,'Enable','off');
    set(H.bandlimits1,'Enable','off');
    set(H.bandlimits2,'Enable','off');
end

% 2nd stim
if ~strcmp(XStimParams.stim_type2, 'File') 	% if not from file
    % if changed, re-set bandwidths
    if XStimParams.bandwidth(2) ~= str2num(get(H.bandwidth2,'String'))
        XStimParams.bandwidth(2) = str2num(get(H.bandwidth2,'String'));
    end
    % if changed, re-set bandlimits
    if (XStimParams.bandlimits(3) ~= str2num(get(H.bandlimits3,'String')) | ...
            XStimParams.bandlimits(4) ~= str2num(get(H.bandlimits4,'String')));  
        XStimParams.bandlimits(3) = str2num(get(H.bandlimits3,'String'));
        XStimParams.bandlimits(4) = str2num(get(H.bandlimits4,'String'));
    end
    set(H.bandwidth2,'Enable','on');
    set(H.bandlimits3,'Enable','on');
    set(H.bandlimits4,'Enable','on');
else
    set(H.bandwidth2,'Enable','off');
    set(H.bandlimits3,'Enable','off');
    set(H.bandlimits4,'Enable','off');
end

% 3rd stim
if ~strcmp(XStimParams.stim_type3, 'File') 	% if not from file
    % if changed, re-set bandwidths
    if XStimParams.bandwidth(3) ~= str2num(get(H.bandwidth3,'String'))
        XStimParams.bandwidth(3) = str2num(get(H.bandwidth3,'String'));
    end
    % if changed, re-set bandlimits
    if (XStimParams.bandlimits(5) ~= str2num(get(H.bandlimits5,'String')) | ...
       XStimParams.bandlimits(6) ~= str2num(get(H.bandlimits6,'String')));  
        XStimParams.bandlimits(5) = str2num(get(H.bandlimits5,'String'));
        XStimParams.bandlimits(6) = str2num(get(H.bandlimits6,'String'));
    end
    set(H.bandwidth3,'Enable','on');
    set(H.bandlimits5,'Enable','on');
    set(H.bandlimits6,'Enable','on');
else
    set(H.bandwidth3,'Enable','off');
    set(H.bandlimits5,'Enable','off');
    set(H.bandlimits6,'Enable','off');
end

% check for visibility of ITD
if XStimParams.ildalone_flag | XStimParams.ABLalone_flag
    set(H.ITD,		'Enable', 'on');
else
    set(H.ITD,		'Enable', 'off');
end

% check for visibility of ILD
if XStimParams.itdalone_flag | XStimParams.ABLalone_flag
    set(H.az,		'Enable', 'on');
    set(H.el,		'Enable', 'on');
else
    set(H.az,		'Enable', 'off');
    set(H.el,		'Enable', 'off');
end

% set ITD, and ABL
XStimParams.curr_ITD = str2num(get(H.ITD,'string'));
XStimParams.curr_ABL = str2num(get(H.ABL,'string'));
XStimParams.test_ISI = str2num(get(H.ISI,'string'));

%check visibility of record datafile
if get(H.recorddata,'Value')
    set(H.recorddata_FN,'enable','on');
else
    set(H.recorddata_FN,'Enable','off');
end

% check out modulation parameters for stim1
XStimParams.mod_depth(1) = str2num(get(H.Twosrc_mod_depth,'String'));
mod_type = get(H.Twosrc_mod_type,'String');
mod_num = get(H.Twosrc_mod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Tone'
        set(H.Twosrc_mod_txt,'String', 'Mod Freq:');
        set(H.Twosrc_mod_freq,'Visible','on');
        set(H.Twosrc_mod_phase,'Visible','on');
        set(H.Twosrc_mod_txtA,'Visible','on');
        set(H.Twosrc_mod_depth,'Visible','on');
        set(H.Twosrc_mod_txtB,'Visible','on');
        set(H.Twosrc_mod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Twosrc_mod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Twosrc_mod_phase,'String'));
    case 'File'
        set(H.Twosrc_mod_pb,'Visible','on');
        set(H.Twosrc_mod_freq,'Visible','off');
        set(H.Twosrc_mod_txtA,'Visible','on');
        set(H.Twosrc_mod_depth,'Visible','on');
        set(H.Twosrc_mod_phase,'Visible','off');
        set(H.Twosrc_mod_txtB,'Visible','off');
        if get(H.Twosrc_mod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.Twosrc_mod_pb,'Value',0);
        end
        set(H.Twosrc_mod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.Twosrc_mod_txt,'String', 'CutOff Freq:');
        set(H.Twosrc_mod_freq,'Visible','on');
        set(H.Twosrc_mod_txtA,'Visible','on');
        set(H.Twosrc_mod_depth,'Visible','on');
        set(H.Twosrc_mod_phase,'Visible','off');
        set(H.Twosrc_mod_pb,'Visible','off');
        set(H.Twosrc_mod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Twosrc_mod_freq,'String'));
        
    case 'None'
        set(H.Twosrc_mod_txt,'String', 'no mod  ');
        set(H.Twosrc_mod_freq,'Visible','off');
        set(H.Twosrc_mod_phase,'Visible','off');
        set(H.Twosrc_mod_pb,'Visible','off');
        set(H.Twosrc_mod_txtB,'Visible','off');
        set(H.Twosrc_mod_txtA,'Visible','off');
        set(H.Twosrc_mod_depth,'Visible','off');
        
    otherwise
end

% check out modulation parameters for stim2
XStimParams.mod_depth(2) = str2num(get(H.Twosrc_mod_depth2,'String'));
mod_type = get(H.Twosrc_mod_type2,'String');
mod_num = get(H.Twosrc_mod_type2,'Value');
XStimParams.mod_type2 = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type2
    case 'Tone'
        set(H.Twosrc_mod_txt2,'String', 'Mod Freq:');
        set(H.Twosrc_mod_freq2,'Visible','on');
        set(H.Twosrc_mod_pb2,'Visible','off');
        set(H.Twosrc_mod_txtA2,'Visible','on');
        set(H.Twosrc_mod_depth2,'Visible','on');
        set(H.Twosrc_mod_phase2,'Visible','on');
        set(H.Twosrc_mod_txtB2,'Visible','on');
        XStimParams.mod_freq(2) = str2num(get(H.Twosrc_mod_freq2,'String'));
        XStimParams.mod_phase(2) = str2num(get(H.Twosrc_mod_phase2,'String'));
    case 'File'
        set(H.Twosrc_mod_pb2,'Visible','on');
        set(H.Twosrc_mod_freq2,'Visible','off');
        set(H.Twosrc_mod_txtA2,'Visible','on');
        set(H.Twosrc_mod_depth2,'Visible','on');
        set(H.Twosrc_mod_phase2,'Visible','off');
        set(H.Twosrc_mod_txtB2,'Visible','off');
        if get(H.Twosrc_mod_pb2,'Value')
            [mod_from_file2, FN.mod2, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.Twosrc_mod_pb2,'Value',0);
        end
        set(H.Twosrc_mod_txt2,'String', FN.mod2);
        
    case 'LP Noise'
        set(H.Twosrc_mod_txt2,'String', 'CutOff Freq:');
        set(H.Twosrc_mod_freq2,'Visible','on');
        set(H.Twosrc_mod_txtA2,'Visible','on');
        set(H.Twosrc_mod_depth2,'Visible','on');
        set(H.Twosrc_mod_pb2,'Visible','off');
        set(H.Twosrc_mod_phase2,'Visible','off');
        set(H.Twosrc_mod_txtB2,'Visible','off');
        XStimParams.mod_freq(2) = str2num(get(H.Twosrc_mod_freq2,'String'));
        
    case 'None'
        set(H.Twosrc_mod_txt2,'String', 'no mod  ');
        set(H.Twosrc_mod_freq2,'Visible','off');
        set(H.Twosrc_mod_phase2,'Visible','off');
        set(H.Twosrc_mod_pb2,'Visible','off');
        set(H.Twosrc_mod_txtB2,'Visible','off');
        set(H.Twosrc_mod_txtA2,'Visible','off');
        set(H.Twosrc_mod_depth2,'Visible','off');
        
    otherwise
end

% check out modulation parameters for stim3
XStimParams.mod_depth(3) = str2num(get(H.Twosrc_mod_depth3,'String'));
mod_type = get(H.Twosrc_mod_type3,'String');
mod_num = get(H.Twosrc_mod_type3,'Value');
XStimParams.mod_type3 = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type3
    case 'Tone'
        set(H.Twosrc_mod_txt3,'String', 'Mod Freq:');
        set(H.Twosrc_mod_freq3,'Visible','on');
        set(H.Twosrc_mod_txtA3,'Visible','on');
        set(H.Twosrc_mod_depth3,'Visible','on');
        set(H.Twosrc_mod_phase3,'Visible','on');
        set(H.Twosrc_mod_pb3,'Visible','off');
        set(H.Twosrc_mod_txtB3,'Visible','on');
        XStimParams.mod_freq(3) = str2num(get(H.Twosrc_mod_freq3,'String'));
        XStimParams.mod_phase(3) = str2num(get(H.Twosrc_mod_phase3,'String'));
    case 'File'
        set(H.Twosrc_mod_pb3,'Visible','on');
        set(H.Twosrc_mod_freq3,'Visible','off');
        set(H.Twosrc_mod_txtA3,'Visible','on');
        set(H.Twosrc_mod_depth3,'Visible','on');
        set(H.Twosrc_mod_phase3,'Visible','off');
        set(H.Twosrc_mod_txtB3,'Visible','off');
        if get(H.Twosrc_mod_pb3,'Value')
            [mod_from_file3, FN.mod3, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.Twosrc_mod_pb3,'Value',0);
        end
        set(H.Twosrc_mod_txt3,'String', FN.mod2);
        
    case 'LP Noise'
        set(H.Twosrc_mod_txt3,'String', 'CutOff Freq:');
        set(H.Twosrc_mod_freq3,'Visible','on');
        set(H.Twosrc_mod_txtA3,'Visible','on');
        set(H.Twosrc_mod_depth3,'Visible','on');
        set(H.Twosrc_mod_pb3,'Visible','off');
        set(H.Twosrc_mod_phase3,'Visible','off');
        set(H.Twosrc_mod_txtB3,'Visible','off');
        XStimParams.mod_freq(3) = str2num(get(H.Twosrc_mod_freq3,'String'));
        
    case 'None'
        set(H.Twosrc_mod_txt3,'String', 'no mod  ');
        set(H.Twosrc_mod_freq3,'Visible','off');
        set(H.Twosrc_mod_pb3,'Visible','off');
        set(H.Twosrc_mod_phase3,'Visible','off');
        set(H.Twosrc_mod_txtB3,'Visible','off');
        set(H.Twosrc_mod_txtA3,'Visible','off');
        set(H.Twosrc_mod_depth3,'Visible','off');
        
    otherwise
end


% calc Vstrength
if get(H.plotVstr,'value')
    tit = 'Vector Strength Plots';
    prompt = {'Enter datapath', 'Enter filename ','Enter vector of modulation frequencies'};
    DefAns = {FN.data_path,[FN.data '.mat'], '[55 75]'};
    
    Ans = inputdlg(prompt,tit,1,DefAns);
    read_data3_minus_spont(Ans{2},Ans{1}, [], 1,str2num(Ans{3}));
    clear tit prompt DefAns Ans
    set(H.plotVstr,'value',0);
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
update_XStim