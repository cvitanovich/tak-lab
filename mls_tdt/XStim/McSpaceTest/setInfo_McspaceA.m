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
if get(H.spacefilepb,'Value')
    [FN.space_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space_std);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
end

XStimParams.curr_ABL = str2num(get(H.mcSpace_ABL,'string'));

% bandwidth & limits
if get(H.stim_type,'Value') ~= 9			% not from file
    XStimParams.bandwidth = str2num(get(H.mcSpace_bandwidth,'String'));
    XStimParams.bandlimits(1) = str2num(get(H.mcSpace_bandlimits1,'String'));
    XStimParams.bandlimits(2) = str2num(get(H.mcSpace_bandlimits2,'String'));
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
set(H.spacedata_FN,'String', FN.data);

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
        str = ['With tone or gammatone stimuli    '; ...
                'b/wn epoch clipping occurs unless '; ...
                'ramppts is non-zero. Larger swings'; ...
                'in seqSPL require higher ramppts. '];
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

% check out modulation parameters for stim
XStimParams.mod_depth(1) = str2num(get(H.Spacetestmod_depth,'String'));
mod_type = get(H.Spacetestmod_type,'String');
mod_num = get(H.Spacetestmod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

set(H.Spacetestmod_txtC, 'Visible','off');
set(H.mcSpace_analVstr,'visible','off');

switch XStimParams.mod_type
    case 'Sq wave'
        set(H.Spacetestmod_txt,'String', 'Freq:');
        set(H.Spacetestmod_freq,'Visible','on');
        set(H.Spacetestmod_phase,'Visible','on');
        set(H.Spacetestmod_txtA,'Visible','on');
        set(H.Spacetestmod_depth,'Visible','on');
        set(H.Spacetestmod_txtB,'Visible','on');
        set(H.Spacetestmod_pb,'Visible','off');
        if stim_val ~= 10   set(H.Spacetestmod_txtC, 'Visible','on');   end
        XStimParams.mod_freq(1) = str2num(get(H.Spacetestmod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Spacetestmod_phase,'String'));
        set(H.mcSpace_analVstr,'visible','on');
    case 'Tone'
        set(H.Spacetestmod_txt,'String', 'Freq:');
        set(H.Spacetestmod_freq,'Visible','on');
        set(H.Spacetestmod_phase,'Visible','on');
        set(H.Spacetestmod_txtA,'Visible','on');
        set(H.Spacetestmod_depth,'Visible','on');
        set(H.Spacetestmod_txtB,'Visible','on');
        set(H.Spacetestmod_pb,'Visible','off');
        if stim_val ~= 10   set(H.Spacetestmod_txtC, 'Visible','on');   end
        XStimParams.mod_freq(1) = str2num(get(H.Spacetestmod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Spacetestmod_phase,'String'));
        set(H.mcSpace_analVstr,'visible','on');
    case 'File'
        set(H.Spacetestmod_pb,'Visible','on');
        set(H.Spacetestmod_freq,'Visible','off');
        set(H.Spacetestmod_txtA,'Visible','on');
        set(H.Spacetestmod_depth,'Visible','on');
        set(H.Spacetestmod_phase,'Visible','off');
        set(H.Spacetestmod_txtB,'Visible','off');
        if stim_val ~= 10   set(H.Spacetestmod_txtC, 'Visible','on');   end
        if get(H.Spacetestmod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.Spacetestmod_pb,'Value',0);
        end
        set(H.Spacetestmod_txt,'String', FN.mod);
    case 'LP Noise'
        set(H.Spacetestmod_txt,'String', 'CutOff Freq:');
        set(H.Spacetestmod_freq,'Visible','on');
        set(H.Spacetestmod_txtA,'Visible','on');
        set(H.Spacetestmod_depth,'Visible','on');
        set(H.Spacetestmod_phase,'Visible','off');
        set(H.Spacetestmod_pb,'Visible','off');
        set(H.Spacetestmod_txtB,'Visible','off');
        if stim_val ~= 10   set(H.Spacetestmod_txtC, 'Visible','on');   end
        XStimParams.mod_freq(1) = str2num(get(H.Spacetestmod_freq,'String'));
    case 'None'
        set(H.Spacetestmod_txt,'String', 'no mod  ');
        set(H.Spacetestmod_freq,'Visible','off');
        set(H.Spacetestmod_phase,'Visible','off');
        set(H.Spacetestmod_pb,'Visible','off');
        set(H.Spacetestmod_txtB,'Visible','off');
        set(H.Spacetestmod_txtA,'Visible','off');
        set(H.Spacetestmod_depth,'Visible','off');
    otherwise
end

%%% values for FOCAL SPLs
XStimParams.epoch_duration(1) = str2num(get(H.epoch_duration,'string'));
XStimParams.focalProb = min([1 str2num(get(H.focalProb,'string'))]);

% reset slider values
for islider = 1:2
    M(islider) = round(get(H.mcSpace_SPLslider(islider),'value'));
end
XStimParams.focalSPL = [min1(M) max1(M)];;
set(H.mcSpace_SPLslider(1),'value',XStimParams.focalSPL(1));
set(H.mcSpace_SPLslider(2),'value',XStimParams.focalSPL(2));
set(H.mcSpace_SPLslider_txt(1),'string',num2str(XStimParams.focalSPL(1)));
set(H.mcSpace_SPLslider_txt(2),'string',num2str(XStimParams.focalSPL(2)));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim