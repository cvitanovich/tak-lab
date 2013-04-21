% SetInfo_Mc2Source
% combines all the other 'set...' scripts into one

%check ILDAloneFlag
% XStimParams.space_flag == 1 for ILDalone
%						 == 0 for fully cued

if get(H.Mc2Src_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.Mc2Src_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.Mc2Src_DynamicRangepb,'BackgroundColor','red');
    set(H.Mc2Src_DynamicRangeTxt,'visible','on');
else
    XStimParams.HiDynamicRange = 0;
    set(H.Mc2Src_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.Mc2Src_DynamicRangepb,'BackgroundColor','yellow');
    set(H.Mc2Src_DynamicRangeTxt,'visible','off');
end

if get(H.Mc2Src_txt(19),'value') % swap to constant S2 and variable S1 probes
    XStimParams.test_type = 'Mc2SourceReversed';
    flag_reverse = 1;
    set(H.Mc2Src_txt(18),'string','S2 constant.    Set S1 (probe) focal-range SPL:');
    set(H.Mc2Src_txt(19),...    
    'BackgroundColor','red',...
    'ForegroundColor','blue',...
    'String', 'norm S1/S2');
    set(H.Mc2Src_txt(6),'string','S1 non-focal SPLs');
    set(H.Mc2Src_txt(16),'String', 'source1 FOCAL total PROB');
    set(H.Mc2Src_txt(17),'String', 'source2 SPL (0 to 30dB)');
    
    set(H.Mc2Src_txt(20), ...
   'BackgroundColor','yellow',...
   'ForegroundColor','blue',...
   'String', 'SetUp SPLs   SWAPPED TO S2 constant, variable S1 probe-level');

else                            % swap to variable S2 and constant S1 probes
    XStimParams.test_type = 'Mc2Source';
    flag_reverse = 0;
    set(H.Mc2Src_txt(18),'string','S1 constant (probe).    Set S2 focal-range SPL:');
    set(H.Mc2Src_txt(19),...    
    'BackgroundColor','yellow',...
    'ForegroundColor','blue',...
    'String', 'swap S1/S2');
    set(H.Mc2Src_txt(6),'string','S2 non-focal SPLs');
    set(H.Mc2Src_txt(16),'String', 'source2 FOCAL total PROB');
    set(H.Mc2Src_txt(17),'String', 'source1 SPL (0 to 30dB)');
    set(H.Mc2Src_txt(20), ...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'String', 'SetUp SPLs');
end
    
if get(H.spaceflag,'value')
    XStimParams.space_flag = 1;
    set(H.spaceflag,'string','push for ABLequal FC')
    set(H.spaceflag,'BackgroundColor','yellow')
    set(H.spacefilepb,'string','Select FC Filter File');   
    set(H.spacefilehdr,'String', 'Fully-cued HRTF File:');
    set(H.spacefile,'string',FN.space_eq);
    XStimParams.ABLequal_flag = 0;
    % check for FC file
    if get(H.spacefilepb,'Value')
        [FN.space_eq,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File *.eq preferred');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',FN.space_eq);
        end
        set(H.spacefilepb,'Value',0);
        FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);
    end
    
else
    XStimParams.space_flag = 0;
    set(H.spaceflag,'string','push for Fully-cued')
    set(H.spaceflag,'BackgroundColor','red')
    set(H.spacefilepb,'String','Select ABLequal FC Filter File');   
    set(H.spacefilehdr,'String', 'ABLequal FC HRTF File:');
    set(H.spacefile,'string',FN.ablequal_eq);
    XStimParams.ABLequal_flag = 1;
    
    % check for ABLequal file
    if get(H.spacefilepb,'Value')
        [FN.ablequal_eq,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File *.eq preferred');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',FN.ablequal_eq);
        end
        set(H.spacefilepb,'Value',0);
        FN.HRTFfiletype(7,1) = testHRTFfiletype(FN.space_path, FN.ablequal_eq);
    end
end

% set numreps
XStimParams.numreps = str2num(get(H.numreps,'String'));

%set trailing and leading silence periods
XStimParams.silence_lead = str2num(get(H.Mc2Source_silence_lead,'string'));
XStimParams.silence_trail = str2num(get(H.Mc2Source_silence_trail,'string'));

% set offsets
XStimParams.offset_el(1) = 0;
XStimParams.offset_az(1) = 0;
XStimParams.offset_el(2) = str2num(get(H.offset_el,'String'));
XStimParams.offset_az(2) = str2num(get(H.offset_az,'String'));

% set multiplication factor for each source
XStimParams.factor(1) = str2num(get(H.Mc2Src_factor,'String'));
XStimParams.factor(2) = str2num(get(H.Mc2Src_factor2,'String'));

% set ramping
XStimParams.ramp = get(H.source1_ramp,'value');

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
if(exist1('H.Mc2Sourcefig') & get(H.inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.inc_testnum,'Value',0);
end
set(H.recorddata_FN,'String',FN.data);

% update stim_type1
stim_type = get(H.Mc2Src_stim_type,'String');
stim_val = get(H.Mc2Src_stim_type,'Value');

if ~strcmp(XStimParams.stim_type, deblank(stim_type(stim_val,:))) | ...
        get(H.Mc2Src_stim_filenamepb,'Value')		% stim_type has changed or pb pushed
    XStimParams.stim_type = deblank(stim_type(stim_val,:));
    set(H.stim_type,'Value',stim_val);
    if(strcmp('File',deblank(XStimParams.stim_type)))
        [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
        set(H.Mc2Src_stim_filename,'String', FN.stim);
        set(H.stim_filename,'String', FN.stim);
    end
    %setinfo_Xstim;
    set(H.Mc2Src_stim_filenamepb,'Value',0);
end

%control visibility of FN vs freq
set(H.Mc2Src_stim_filenamepb,'Visible','on');
set(H.Mc2Src_stim_filename,'Visible','on');
set(H.Mc2Src_freq_txt,'Visible','off');
set(H.Mc2Src_freq,'Visible','off');    

switch deblank(XStimParams.stim_type)
    case 'File'
        set(H.Mc2Src_stim_filename,'Enable','on');
        set(H.stim_filename,'Enable','on');
    case {'Gammatone', 'Tone'}
        set(H.Mc2Src_freq_txt,'Visible','on');
        set(H.Mc2Src_freq,'Visible','on');    
        set(H.Mc2Src_stim_filenamepb,'Visible','off');
        set(H.Mc2Src_stim_filename,'Visible','off');
    otherwise
        set(H.Mc2Src_stim_filename,'Enable','off');
        set(H.stim_filename,'Enable','off');
end

% update stim_type2
stim_type = get(H.Mc2Src_stim_type2,'String');
stim_val = get(H.Mc2Src_stim_type2,'Value');

if ~strcmp(XStimParams.stim_type2,deblank(stim_type(stim_val,:))) | ...
        get(H.Mc2Src_stim_filename2pb,'Value')		% stim_type2 has changed
    XStimParams.stim_type2 = deblank(stim_type(stim_val,:));
    set(H.stim_type2,'Value',stim_val);
    if(strcmp('File',deblank(XStimParams.stim_type2)))
        [stim_from_file2, FN.stim2, FN.stim_path2] = Get_Stimulus_File(XStimParams, FN.stim_path2);
        set(H.Mc2Src_stim_filename2,'String', FN.stim2);
        set(H.stim_filename2,'String', FN.stim2);
    end
    set(H.Mc2Src_stim_filename2pb,'Value',0);
end
%control visibility of FN vs freq
set(H.Mc2Src_stim_filename2pb,'Visible','on');
set(H.Mc2Src_stim_filename2,'Visible','on');
set(H.Mc2Src_freq2_txt,'Visible','off');
set(H.Mc2Src_freq2,'Visible','off');    

switch deblank(XStimParams.stim_type2)
    case 'File'
        set(H.Mc2Src_stim_filename2,'Enable','on');
        set(H.stim_filename2,'Enable','on');
    case {'Gammatone', 'Tone'}
        set(H.Mc2Src_freq2_txt,'Visible','on');
        set(H.Mc2Src_freq2,'Visible','on');    
        set(H.Mc2Src_stim_filename2pb,'Visible','off');
        set(H.Mc2Src_stim_filename2,'Visible','off');
    otherwise
        set(H.Mc2Src_stim_filename2,'Enable','off');
        set(H.stim_filename2,'Enable','off');
end

% set dur based on larger of stim1, stim2 
set(H.DUR,'Enable','on');
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
temp = zeros(4,1);
temp(1) = XStimParams.curr_stimdur;
if (get(H.Mc2Src_stim_type,'Value') == 9)       % from file
    set(H.DUR,'Enable','off');
    clear dir
    D = dir([FN.stim_path FN.stim]);
    temp(2) = round(1000 * (D.bytes/4) / TDT.Fs);
end
if (get(H.Mc2Src_stim_type2,'Value') == 9)       % from file
    set(H.DUR,'Enable','off');
    clear dir
    D = dir([FN.stim_path2 FN.stim2]);
    temp(3) = round(1000 * (D.bytes/4) / TDT.Fs);
end
XStimParams.curr_stimdur = max1(temp);

% repeat mode params
%%% values for FOCAL SPLs
XStimParams.epoch_duration(1) = str2num(get(H.epoch_duration,'string'));
XStimParams.focalProb = min([1 str2num(get(H.focalProb,'string'))]);
XStimParams.source1_SPL = str2num(get(H.source1_SPL,'String'));
XStimParams.source1_ISI = str2num(get(H.source1_ISI,'String'));

% reset slider values
for islider = 1:2
    M(islider) = round(get(H.Mc2Source_SPLslider(islider),'value'));
end
XStimParams.focalSPL = [min1(M) max1(M)];;
set(H.Mc2Source_SPLslider(1),'value',XStimParams.focalSPL(1));
set(H.Mc2Source_SPLslider(2),'value',XStimParams.focalSPL(2));
set(H.Mc2Source_SPLslider_txt(1),'string',num2str(XStimParams.focalSPL(1)));
set(H.Mc2Source_SPLslider_txt(2),'string',num2str(XStimParams.focalSPL(2)));
if ~diff(XStimParams.focalSPL)
    XStimParams.focalProb =1;
end
set(H.focalProb,'string',num2str(XStimParams.focalProb));

if get(H.Mc2Src_repeatmodepb,'value')
    XStimParams.repeatmode = 1;
    set(H.Mc2Src_repeatmodepb,'string','push for continuous mode')
    set(H.Mc2Src_repeatmodepb,'BackgroundColor','red')
    set(H.Mc2Src_onEpochs,'visible','on')
    XStimParams.onEpochs = str2num(get(H.Mc2Src_onEpochs,'string'));
    
    set(H.Mc2Src_offEpochs,'visible','on')
    XStimParams.offEpochs = str2num(get(H.Mc2Src_offEpochs,'string'));
    
    set(H.McSrc_probeFlag(1),'visible','on');
    set(H.McSrc_probeFlag(2),'visible','on');
    XStimParams.probeFlag = [get(H.McSrc_probeFlag(1),'value') get(H.McSrc_probeFlag(2),'value')];
    %XStimParams.probeFlag = [1 get(H.McSrc_probeFlag(2),'value')];
    if XStimParams.source1_ISI < (XStimParams.onEpochs + XStimParams.offEpochs);
        XStimParams.offEpochs = XStimParams.offEpochs - mod((XStimParams.onEpochs+XStimParams.offEpochs),XStimParams.source1_ISI);
    else
        XStimParams.offEpochs = XStimParams.source1_ISI - XStimParams.onEpochs;  
    end
    set(H.Mc2Src_offEpochs,'string',num2str(XStimParams.offEpochs));
    
    for i = 9:15
        set(H.Mc2Src_txt(i),'visible','on')
    end
    
    temp = (XStimParams.onEpochs + XStimParams.offEpochs) * XStimParams.epoch_duration(1);
    temp2 = floor((XStimParams.curr_stimdur- XStimParams.epoch_duration(1))/temp);
    XStimParams.curr_stimdur = max1([temp floor(XStimParams.curr_stimdur/temp)*temp])+ XStimParams.epoch_duration(1);
else
    XStimParams.repeatmode = 0;
    set(H.Mc2Src_repeatmodepb,'string','push for repeat mode')
    set(H.Mc2Src_repeatmodepb,'BackgroundColor','yellow')
    set(H.Mc2Src_onEpochs,'visible','off')
    set(H.Mc2Src_offEpochs,'visible','off')
    for i = [9:13 15]
        set(H.Mc2Src_txt(i),'visible','off')
    end    
    set(H.McSrc_probeFlag(2),'visible','off');
    XStimParams.probeFlag(2) = 0;
    set(H.McSrc_probeFlag(2),'value',0);
    XStimParams.curr_stimdur = str2num(get(H.DUR,'string'));
end

XStimParams.curr_stimdur = min1([XStimParams.curr_stimdur 30000]);
set(H.DUR,'string',num2str(XStimParams.curr_stimdur));

% get SPL increment
XStimParams.SPLinc(1) = str2num(get(H.Mc2Source_SPLinc, 'String'));

% bandwidths and bandlimits
% 1st stim
if ismember(get(H.Mc2Src_stim_type,'Value'), [5 6 7 9])	% if not from file
    % if changed, re-set bandlimits
    if (XStimParams.bandlimits(1) ~= str2num(get(H.bandlimits1,'String')) | ...
            XStimParams.bandlimits(2) ~= str2num(get(H.bandlimits2,'String')));  
        XStimParams.bandlimits(1) = str2num(get(H.bandlimits1,'String'));
        XStimParams.bandlimits(2) = str2num(get(H.bandlimits2,'String'));
    end
    set(H.bandlimits1,'Visible','on');
    set(H.bandlimits2,'Visible','on');
    set(H.bandlimits_txt1,'Visible','on');
else
    set(H.bandlimits_txt1,'Visible','off');
    set(H.bandlimits1,'Visible','off');
    set(H.bandlimits2,'Visible','off');
end

% 2nd stim
if ismember(get(H.Mc2Src_stim_type2,'Value'), [5 6 7 9])	% if not from file
    % if changed, re-set bandlimits
    if (XStimParams.bandlimits(3) ~= str2num(get(H.bandlimits3,'String')) | ...
            XStimParams.bandlimits(4) ~= str2num(get(H.bandlimits4,'String')));  
        XStimParams.bandlimits(3) = str2num(get(H.bandlimits3,'String'));
        XStimParams.bandlimits(4) = str2num(get(H.bandlimits4,'String'));
    end
    set(H.bandlimits_txt2,'Visible','on');
    set(H.bandlimits3,'Visible','on');
    set(H.bandlimits4,'Visible','on');
else
    set(H.bandlimits_txt2,'Visible','off');
    set(H.bandlimits3,'Visible','off');
    set(H.bandlimits4,'Visible','off');
end


% set ABL
XStimParams.curr_ABL = round(str2num(get(H.ABL,'string')));
set(H.ABL,'string',num2str(XStimParams.curr_ABL));
XStimParams.test_ISI = str2num(get(H.ISI,'string'));

%check visibility of record datafile
if get(H.recorddata,'Value')
    set(H.recorddata_FN,'enable','on');
else
    set(H.recorddata_FN,'Enable','off');
end

% check out modulation parameters for stim1
XStimParams.mod_depth(1) = str2num(get(H.Mc2Src_mod_depth,'String'));
mod_type = get(H.Mc2Src_mod_type,'String');
mod_num = get(H.Mc2Src_mod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Tone'
        set(H.Mc2Src_mod_txt,'String', 'Mod Freq:');
        set(H.Mc2Src_mod_freq,'Visible','on');
        set(H.Mc2Src_mod_phase,'Visible','on');
        set(H.Mc2Src_mod_txtA,'Visible','on');
        set(H.Mc2Src_mod_depth,'Visible','on');
        set(H.Mc2Src_mod_txtB,'Visible','on');
        set(H.Mc2Src_mod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Mc2Src_mod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Mc2Src_mod_phase,'String'));
    case 'File'
        set(H.Mc2Src_mod_pb,'Visible','on');
        set(H.Mc2Src_mod_freq,'Visible','off');
        set(H.Mc2Src_mod_txtA,'Visible','on');
        set(H.Mc2Src_mod_depth,'Visible','on');
        set(H.Mc2Src_mod_phase,'Visible','off');
        set(H.Mc2Src_mod_txtB,'Visible','off');
        if get(H.Mc2Src_mod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.Mc2Src_mod_pb,'Value',0);
        end
        set(H.Mc2Src_mod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.Mc2Src_mod_txt,'String', 'CutOff Freq:');
        set(H.Mc2Src_mod_freq,'Visible','on');
        set(H.Mc2Src_mod_txtA,'Visible','on');
        set(H.Mc2Src_mod_depth,'Visible','on');
        set(H.Mc2Src_mod_phase,'Visible','off');
        set(H.Mc2Src_mod_pb,'Visible','off');
        set(H.Mc2Src_mod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Mc2Src_mod_freq,'String'));
        
    case 'None'
        set(H.Mc2Src_mod_txt,'String', 'no mod  ');
        set(H.Mc2Src_mod_freq,'Visible','off');
        set(H.Mc2Src_mod_phase,'Visible','off');
        set(H.Mc2Src_mod_pb,'Visible','off');
        set(H.Mc2Src_mod_txtB,'Visible','off');
        set(H.Mc2Src_mod_txtA,'Visible','off');
        set(H.Mc2Src_mod_depth,'Visible','off');
        
    otherwise
end

% check out modulation parameters for stim2
XStimParams.mod_depth(2) = str2num(get(H.Mc2Src_mod_depth2,'String'));
mod_type = get(H.Mc2Src_mod_type2,'String');
mod_num = get(H.Mc2Src_mod_type2,'Value');
XStimParams.mod_type2 = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type2
    case 'Tone'
        set(H.Mc2Src_mod_txt2,'String', 'Mod Freq:');
        set(H.Mc2Src_mod_freq2,'Visible','on');
        set(H.Mc2Src_mod_pb2,'Visible','off');
        set(H.Mc2Src_mod_txtA2,'Visible','on');
        set(H.Mc2Src_mod_depth2,'Visible','on');
        set(H.Mc2Src_mod_phase2,'Visible','on');
        set(H.Mc2Src_mod_txtB2,'Visible','on');
        XStimParams.mod_freq(2) = str2num(get(H.Mc2Src_mod_freq2,'String'));
        XStimParams.mod_phase(2) = str2num(get(H.Mc2Src_mod_phase2,'String'));
    case 'File'
        set(H.Mc2Src_mod_pb2,'Visible','on');
        set(H.Mc2Src_mod_freq2,'Visible','off');
        set(H.Mc2Src_mod_txtA2,'Visible','on');
        set(H.Mc2Src_mod_depth2,'Visible','on');
        set(H.Mc2Src_mod_phase2,'Visible','off');
        set(H.Mc2Src_mod_txtB2,'Visible','off');
        if get(H.Mc2Src_mod_pb2,'Value')
            [mod_from_file2, FN.mod2, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.Mc2Src_mod_pb2,'Value',0);
        end
        set(H.Mc2Src_mod_txt2,'String', FN.mod2);
        
    case 'LP Noise'
        set(H.Mc2Src_mod_txt2,'String', 'CutOff Freq:');
        set(H.Mc2Src_mod_freq2,'Visible','on');
        set(H.Mc2Src_mod_txtA2,'Visible','on');
        set(H.Mc2Src_mod_depth2,'Visible','on');
        set(H.Mc2Src_mod_pb2,'Visible','off');
        set(H.Mc2Src_mod_phase2,'Visible','off');
        set(H.Mc2Src_mod_txtB2,'Visible','off');
        XStimParams.mod_freq(2) = str2num(get(H.Mc2Src_mod_freq2,'String'));
        
    case 'None'
        set(H.Mc2Src_mod_txt2,'String', 'no mod  ');
        set(H.Mc2Src_mod_freq2,'Visible','off');
        set(H.Mc2Src_mod_phase2,'Visible','off');
        set(H.Mc2Src_mod_pb2,'Visible','off');
        set(H.Mc2Src_mod_txtB2,'Visible','off');
        set(H.Mc2Src_mod_txtA2,'Visible','off');
        set(H.Mc2Src_mod_depth2,'Visible','off');
        
    otherwise
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
update_XStim