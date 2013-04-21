% SetInfo_PE_env

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
XStimParams.silence_lead = str2num(get(H.PE_env_silence_lead,'string'));
XStimParams.silence_trail = str2num(get(H.PE_env_silence_trail,'string'));

% set offsets
XStimParams.offset_el(1) = 0;
XStimParams.offset_az(1) = 0;
XStimParams.offset_el(2) = str2num(get(H.offset_el,'String'));
XStimParams.offset_az(2) = str2num(get(H.offset_az,'String'));

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
if(exist1('H.PE_envfig') & get(H.inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.inc_testnum,'Value',0);
end
set(H.recorddata_FN,'String',FN.data);


% set dur 
set(H.DUR,'Enable','on');
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));

% set ABL
XStimParams.curr_ABL = round(str2num(get(H.ABL,'string')));
set(H.ABL,'string',num2str(XStimParams.curr_ABL));

% set ISI
%XStimParams.test_ISI = str2num(get(H.ISI,'string'));

% set delay(s)
XStimParams.delay = round(str2num(get(H.PE_env_delay,'string')));
set(H.PE_env_delay,'string',num2str(XStimParams.delay));

% set gammatone cF(s)
XStimParams.cF = round(str2num(get(H.PE_env_cF,'string')));
set(H.PE_env_cF,'string',num2str(XStimParams.cF));


%check visibility of record datafile
if get(H.recorddata,'Value')
    set(H.recorddata_FN,'enable','on');
else
    set(H.recorddata_FN,'Enable','off');
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
update_XStim