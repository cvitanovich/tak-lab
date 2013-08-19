% SetInfo_Mcspace

% uses .eq HRIRs, therefore make sure earphone filters are also selected
XStimParams.ephone_flag = 0;
set(H.ephoneuseit,'Value',0);
XStimParams.space_flag = 1;
XStimParams.ildalone_flag = 0;
XStimParams.itdalone_flag = 0;
XStimParams.ABLalone_flag = 0;

if get(H.mcSpace_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.mcSpace_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.mcSpace_DynamicRangepb,'BackgroundColor','red');
    set(H.mcSpace_DynamicRangeTxt,'visible','on');
else
    XStimParams.HiDynamicRange = 0;
    set(H.mcSpace_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.mcSpace_DynamicRangepb,'BackgroundColor','yellow');
    set(H.mcSpace_DynamicRangeTxt,'visible','off');
end

% check for fully-cued pushbutton (file-picker)
if get(H.mcSpace_filepb,'Value')
    [FN.space_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File');
    if(FN.space_path ~= 0)
        set(H.mcSpace_file,'String',FN.space_std);
    end
    set(H.mcSpace_filepb,'Value',0);
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
end

XStimParams.curr_ABL = str2num(get(H.mcSpace_ABL,'string'));

% bandwidth & limits
if get(H.stim_type,'Value') ~= 9			% not from file
    XStimParams.bandwidth = str2num(get(H.mcSpace_bandwidth,'String'));
    XStimParams.bandlimits(1) = str2num(get(H.mcSpace_bandlimits(1),'String'));
    XStimParams.bandlimits(2) = str2num(get(H.mcSpace_bandlimits(2),'String'));
end

% set numreps
XStimParams.numreps = str2num(get(H.mcSpace_numreps,'String'));

% set ramppts
temp = str2num(get(H.mcSpace_ramppts,'string'));
if mod(temp,2)    temp = temp+1;  end 
set(H.mcSpace_ramppts,'string',num2str(temp));
XStimParams.ramppts = temp;
clear temp

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
if(exist1('H.mcSpacefig') & get(H.mcSpace_inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.mcSpace_inc_testnum,'Value',0);
end
set(H.mcSpace_data_FN,'String', FN.data);

% update display of stim filename
stim_type = get(H.mcSpace_stim_type,'String');
stim_val = get(H.mcSpace_stim_type,'Value');
set(H.stim_type,'Value',stim_val);
XStimParams.stim_type = deblank(stim_type(stim_val,:));

clear stim_type

set(H.stim_filename,'Enable','off');
set(H.mcSpace_stimFN,'Enable','off');
set(H.mcSpace_freqtxt,'visible','off');
set(H.mcSpace_freq,'visible','off');

switch stim_val
    case {1 2}
        set(H.mcSpace_freqtxt,'visible','on');
        set(H.mcSpace_freq,'visible','on');
        XStimParams.curr_freq = str2num(get(H.mcSpace_freq,'string'));
        str = ['w/ tone or gammatones   '; ...
                'clips b/wn epoch, unless'; ...
                'ramppts is non-zero.    '; ...
                'high SPL range -> ^^ pts'];
        set(H.mcSpace_rampptstxt(2),'string',str);
        set(H.mcSpace_rampptstxt(2),'ForeGroundColor',[0 .5 .5])
    case 9								% from file
        set(H.stim_filename,'String',FN.stim);
        set(H.mcSpace_stimFN,'String',FN.stim);
        set(H.mcSpace_stimFN,'Enable','on');
        set(H.mcSpace_stimFNpb,'Enable','on');
        set(H.mcSpace_repeatmode,'value',1);
    otherwise
        str = ['ramppts is traditionally set'; ...
                'to zero for this stimtype.  '];
        set(H.mcSpace_rampptstxt(2),'string',str);
        set(H.mcSpace_rampptstxt(2),'ForeGroundColor','blue')
end

if get(H.mcSpace_stimFNpb,'value')
    [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
    set(H.mcSpace_stimFN,'Enable','on');
    D = dir([FN.stim_path FN.stim]);
    temp = (D.bytes/4) - (XStimParams.silence_lead + XStimParams.silence_trail)*(TDT.Fs/1000) - 254;
    XStimParams.curr_stimdur = round(1000 * temp / TDT.Fs);
    set(H.mcSpace_stimFN,'String',FN.stim);
end

XStimParams.repeatmode = get(H.mcSpace_repeatmode,'value');

% update stimulus duration
if stim_val == 9
    set(H.mcSpace_DUR, 'String',num2str(XStimParams.curr_stimdur));
else
    XStimParams.curr_stimdur = str2num(get(H.mcSpace_DUR, 'String'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempstr = ['SPLtest  ';'AMdepth  '; 'AMdepthdB'; 'AMperiod ';...
        'AMperiodB'; 'ILDtest  '; 'ILDtest L'; 'ILDtest R'];
% is this a new type of test?"
newtype = 0;
for iType = 1:8
    if get(H.mcSpace_McTestType(iType),'value')
        if ~strcmp(XStimParams.McTestType,tempstr(iType,:))
            newtype=iType;
        end
    end
end
if newtype
    for iType = setdiff(1:8,newtype)
        set(H.mcSpace_McTestType(iType),'value',0);
    end
end

% get McTestType
McTestType = 'dunno ';
for iType = 1:8
    if get(H.mcSpace_McTestType(iType),'value')
        McTestType = tempstr(iType,:);
    end
end

if ~strcmp(McTestType,'ILDtest R') & ~strcmp(McTestType,'ILDtest L')
    set(H.mcSpace_TestTypetxt(4), 'visible','off'); 
    set(H.mcSpace_offEpochs,'visible','off');
end

clear M
switch McTestType
    case 'SPLtest  '
        set(H.mcSpace_DURtxt(2),'String', '32K max');
        if ~strcmp(XStimParams.McTestType,McTestType)
            set(H.mcSpace_mod_txt(4),'visible','on');          % mod phase text
            set(H.mcSpace_mod_phase,'visible','on');           % mod phase
            set(H.mcSpace_mod_txt(1),'visible','on');          % mod freq text
            set(H.mcSpace_mod_freq ,'visible','on');           % mod freq
            set(H.mcSpace_mod_txt(2),'visible','on');          % mod depth text
            set(H.mcSpace_mod_depth ,'visible','on');          % mod depth
            set(H.mcSpace_epoch_durationtxt,'visible','on');   % epoch_dur text
            set(H.mcSpace_epoch_duration,'visible','on');      % epoch_dur
            
            set(H.mcSpace_SPLtxt(3),'visible','off');
            set(H.mcSpace_SPLtxt(1),'string','SetUp SPLs');
            set(H.mcSpace_SPLtxt(2),'string','focal-range SPL');
            % reset slider values
            M = [0 30];
            for islider = 1:2
                set(H.mcSpace_SPLslider(islider),'Min',M(1),'Max',M(2),'SliderStep',[1/31 5/31]);
                set(H.mcSpace_SPLslider(islider),'value',M(islider));
                set(H.mcSpace_SPLslider_txt(islider),'string',num2str(M(islider)));
            end
        else
            % get slider values
            for islider = 1:2
                M(islider) = round(get(H.mcSpace_SPLslider(islider),'value'));
            end
            XStimParams.focalSPL = [min1(M) max1(M)];
            for islider = 1:2
                set(H.mcSpace_SPLslider(islider),'value',XStimParams.focalSPL(islider));
                set(H.mcSpace_SPLslider_txt(islider),'string',num2str(XStimParams.focalSPL(islider)));
            end
        end
    case {'AMdepth  ', 'AMdepthdB'}
        set(H.mcSpace_DURtxt(2),'String', '32K max');
        if ~strcmp(XStimParams.McTestType,McTestType)
            set(H.mcSpace_mod_txt(4),'visible','off');          % mod phase text
            set(H.mcSpace_mod_phase,'visible','off');           % mod phase
            XStimParams.mod_phase(1) = pi/2;                    % mod phase
            set(H.mcSpace_mod_txt(1),'visible','on');           % mod freq text
            set(H.mcSpace_mod_freq ,'visible','on');            % mod freq
            set(H.mcSpace_mod_txt(2),'visible','off');          % mod depth text
            set(H.mcSpace_mod_depth ,'visible','off');          % mod depth
            set(H.mcSpace_epoch_durationtxt,'visible','off');   % epoch_dur text
            set(H.mcSpace_epoch_duration,'visible','off');      % epoch_dur
            
            set(H.mcSpace_SPLtxt(3),'visible','off');
            set(H.mcSpace_SPLtxt(1),'string','SetUp AM depths');
            set(H.mcSpace_SPLtxt(2),'string','focal-range AMdepth');
            
            M = [.1 1];
            for islider = 1:2
                set(H.mcSpace_SPLslider(islider),'Min',M(1),'Max',M(2),'SliderStep',[1/9 1/9]);
                set(H.mcSpace_SPLslider(islider),'value',M(islider));
                set(H.mcSpace_SPLslider_txt(islider),'string',num2str(M(islider)));
            end
        else
            % reset slider values to handle AMdepths - retain values in
            % XStimParams as if SPLs
            for islider = 1:2
                M(islider) = get(H.mcSpace_SPLslider(islider),'value');
            end
            XStimParams.focalSPL = [min1(M) max1(M)];
            for islider = 1:2
                set(H.mcSpace_SPLslider(islider),'value',XStimParams.focalSPL(islider));
                set(H.mcSpace_SPLslider_txt(islider),'string',num2str(XStimParams.focalSPL(islider)));
            end
        end
    case {'AMperiod ', 'AMperiodB'}
        set(H.mcSpace_DURtxt(2),'String', '32K max');
        if ~strcmp(XStimParams.McTestType,McTestType)
            set(H.mcSpace_mod_txt(4),'visible','off');          % mod phase text
            set(H.mcSpace_mod_phase,'visible','off');           % mod phase
            XStimParams.mod_phase(1) = pi/2;                    % mod phase
            set(H.mcSpace_mod_txt(1),'visible','off');          % mod freq text
            set(H.mcSpace_mod_freq ,'visible','off');           % mod freq
            set(H.mcSpace_mod_txt(2),'visible','on');           % mod depth text
            set(H.mcSpace_mod_depth ,'visible','on');           % mod depth
            set(H.mcSpace_epoch_durationtxt,'visible','off');   % epoch_dur text
            set(H.mcSpace_epoch_duration,'visible','off');      % epoch_dur
            
            set(H.mcSpace_SPLtxt(1),'string','SetUp AM periods');
            set(H.mcSpace_SPLtxt(2),'string','focal-range AMperiod');
            set(H.mcSpace_SPLtxt(3),'visible','on');
            
            Mperiods = [150 158 167 176 187 200 214 231 250 273 300 333 375 429 500 600 750 1000 1500];
            Mfreqs = round((1000./[Mperiods(1) Mperiods(end)]./(1/30)));
            M = [1 length(Mperiods)];
            for islider = 1:2
                set(H.mcSpace_SPLslider(islider),'Min',M(1),'Max',M(2),'SliderStep',[1 1]/M(2));
                set(H.mcSpace_SPLslider(islider),'value',M(islider));
                set(H.mcSpace_SPLslider_txt(islider),'string',num2str(Mfreqs(islider)));
            end
        else
            % reset slider values to handle AMperiods - retain values in
            % XStimParams as if indices to Mperiods
            Mperiods = [150 158 167 176 187 200 214 231 250 273 300 333 375 429 500 600 750 1000 1500];
            Mfreqs = round((1000./Mperiods./(1/30)));       % where dt = 1/30;
            for islider = 1:2
                M(islider) = round(get(H.mcSpace_SPLslider(islider),'value'));
            end
            XStimParams.focalSPL = [min1(M) max1(M)];
            for islider = 1:2
                set(H.mcSpace_SPLslider(islider),'value',XStimParams.focalSPL(islider));
                set(H.mcSpace_SPLslider_txt(islider),'string',num2str(Mfreqs(XStimParams.focalSPL(islider))));
            end
        end
        
    case {'ILDtest  ', 'ILDtest R', 'ILDtest L'}
        set(H.mcSpace_DURtxt(2),'String', '24K max');
        if ~strcmp(XStimParams.McTestType,McTestType)
            set(H.mcSpace_mod_txt(4),'visible','off');          % mod phase text
            set(H.mcSpace_mod_phase,'visible','off');           % mod phase
            set(H.mcSpace_mod_txt(1),'visible','off');          % mod freq text
            set(H.mcSpace_mod_freq ,'visible','off');           % mod freq
            set(H.mcSpace_mod_txt(2),'visible','off');          % mod depth text
            set(H.mcSpace_mod_depth ,'visible','off');          % mod depth
            set(H.mcSpace_epoch_durationtxt,'visible','on');   % epoch_dur text
            set(H.mcSpace_epoch_duration,'visible','on');      % epoch_dur
            
            set(H.mcSpace_SPLtxt(3),'visible','off');
            set(H.mcSpace_SPLtxt(1),'string','SetUp ILDs');
            set(H.mcSpace_SPLtxt(2),'string','focal-range ILD');
            % reset slider values
            M = [-20 20];
            for islider = 1:2
                set(H.mcSpace_SPLslider(islider),'Min',M(1),'Max',M(2),'SliderStep',[1/41 5/41]);
                set(H.mcSpace_SPLslider(islider),'value',M(islider));
                set(H.mcSpace_SPLslider_txt(islider),'string',num2str(M(islider)));
            end
        else
            % get slider values
            for islider = 1:2
                M(islider) = round(get(H.mcSpace_SPLslider(islider),'value'));
            end
            XStimParams.focalSPL = [min1(M) max1(M)];
            for islider = 1:2
                set(H.mcSpace_SPLslider(islider),'value',XStimParams.focalSPL(islider));
                set(H.mcSpace_SPLslider_txt(islider),'string',num2str(XStimParams.focalSPL(islider)));
            end
        end
        if strcmp(McTestType,'ILDtest R') | strcmp(McTestType,'ILDtest L')
            set(H.mcSpace_TestTypetxt(4), 'visible','on'); 
            set(H.mcSpace_offEpochs,'visible','on');
            XStimParams.offEpochs = str2num(get(H.mcSpace_offEpochs,'string'));
        end
        
end

XStimParams.McTestType = McTestType;
%%% values for FOCAL SPLs
XStimParams.epoch_duration(1) = str2num(get(H.mcSpace_epoch_duration,'string'));
XStimParams.focalProb = min([1 str2num(get(H.mcSpace_focalProb,'string'))]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check out modulation parameters for stim
XStimParams.mod_depth(1) = str2num(get(H.mcSpace_mod_depth,'String'));
XStimParams.mod_freq(1) = str2num(get(H.mcSpace_mod_freq,'String'));

if strcmp(McTestType,'SPLtest  ')
    XStimParams.mod_phase(1) = str2num(get(H.mcSpace_mod_phase,'String'));
   
    mod_type = get(H.mcSpace_mod_type,'String');
    mod_num = get(H.mcSpace_mod_type,'Value');
    XStimParams.mod_type = deblank(mod_type(mod_num,:));
    
    set(H.mcSpace_mod_txt(3), 'Visible','off');
    set(H.mcSpace_analVstr,'visible','off');
    
    if strcmp( XStimParams.mod_type,'Sq wave')
        set(H.mcSpace_mod_txt(1),'String', 'Freq:');
        set(H.mcSpace_mod_freq,'Visible','on');
        set(H.mcSpace_mod_phase,'Visible','on');
        set(H.mcSpace_mod_txt(2),'Visible','on');
        set(H.mcSpace_mod_depth,'Visible','on');
        set(H.mcSpace_mod_txt(3),'Visible','on');
        set(H.mcSpace_mod_pb,'Visible','off');
        if stim_val ~= 10   set(H.mcSpace_mod_txt(3), 'Visible','on');   end
        set(H.mcSpace_analVstr,'visible','on');
    elseif strcmp( XStimParams.mod_type, 'Tone')
        set(H.mcSpace_mod_txt(1),'String', 'Freq:');
        set(H.mcSpace_mod_freq,'Visible','on');
        set(H.mcSpace_mod_phase,'Visible','on');
        set(H.mcSpace_mod_txt(2),'Visible','on');
        set(H.mcSpace_mod_depth,'Visible','on');
        set(H.mcSpace_mod_txt(3),'Visible','on');
        set(H.mcSpace_mod_pb,'Visible','off');
        if stim_val ~= 10   set(H.mcSpace_mod_txt(3), 'Visible','on');   end
        set(H.mcSpace_analVstr,'visible','on');
    elseif strcmp( XStimParams.mod_type, 'File')
        set(H.mcSpace_mod_txt(1),'String', FN.mod);
        set(H.mcSpace_mod_pb,'Visible','on');
        set(H.mcSpace_mod_freq,'Visible','off');
        set(H.mcSpace_mod_txt(2),'Visible','on');
        set(H.mcSpace_mod_depth,'Visible','on');
        set(H.mcSpace_mod_phase,'Visible','off');
        set(H.mcSpace_mod_txt(3),'Visible','off');
        if stim_val ~= 10   set(H.mcSpace_mod_txt(3), 'Visible','on');   end
        if get(H.mcSpace_mod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.mcSpace_mod_pb,'Value',0);
        end
    elseif strcmp( XStimParams.mod_type, 'LP Noise')
        set(H.mcSpace_mod_txt(1),'String', 'CutOff Freq:');
        set(H.mcSpace_mod_freq,'Visible','on');
        set(H.mcSpace_mod_txt(2),'Visible','on');
        set(H.mcSpace_mod_depth,'Visible','on');
        set(H.mcSpace_mod_phase,'Visible','off');
        set(H.mcSpace_mod_pb,'Visible','off');
        set(H.mcSpace_mod_txt(3),'Visible','off');
        if stim_val ~= 10   set(H.mcSpace_mod_txt(3), 'Visible','on');   end
    elseif strcmp( XStimParams.mod_type, 'None')
        set(H.mcSpace_mod_txt(1),'String', 'no mod  ');
        set(H.mcSpace_mod_freq,'Visible','off');
        set(H.mcSpace_mod_phase,'Visible','off');
        set(H.mcSpace_mod_pb,'Visible','off');
        set(H.mcSpace_mod_txt(2),'Visible','off');
        set(H.mcSpace_mod_txt(3),'Visible','off');
        set(H.mcSpace_mod_depth,'Visible','off');
    end
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim