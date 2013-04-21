% SetInfo_MaskedSpace
% combines all the other 'set...' scripts into one


% uses .eq HRIRs
XStimParams.ephone_flag = 0;
XStimParams.buildOnly = 0;

XStimParams.test_ISI = str2num(get(H.MaskedSpace_ISI,'string'));
set(H.Xstim_test_ISI,'string',num2str(XStimParams.test_ISI));

if get(H.MaskedSpace_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.MaskedSpace_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.MaskedSpace_DynamicRangepb,'BackgroundColor','red');
    set(H.MaskedSpace_DynamicRangeTxt,'visible','on');
else
    XStimParams.HiDynamicRange = 0;
    set(H.MaskedSpace_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.MaskedSpace_DynamicRangepb,'BackgroundColor','yellow');
    set(H.MaskedSpace_DynamicRangeTxt,'visible','off');
end

% check additional dB for masker and probe
XStimParams.mask_add_dB = str2num(get(H.MaskedSpace_mask_add_dB,'string'));
XStimParams.probe_add_dB = str2num(get(H.MaskedSpace_probe_add_dB,'string'));
XStimParams.mask_add_dB = min1([20 XStimParams.mask_add_dB]);
XStimParams.probe_add_dB = min1([20 XStimParams.probe_add_dB]);
set(H.MaskedSpace_mask_add_dB,'string',num2str(XStimParams.mask_add_dB))
set(H.MaskedSpace_probe_add_dB,'string',num2str(XStimParams.probe_add_dB))

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

XStimParams.curr_ABL = str2num(get(H.MaskedSpace_ABL,'string'));

% set numreps
XStimParams.numreps = str2num(get(H.MaskedSpace_numreps,'String'));

%Set space PickerFlag - once engaged, pick probe location first, then picker clears and
%allows picking of masker locas
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
set(H.pickerfig,'Name','Probe Space Picker')

% increment test number
if(exist1('H.MaskedSpacefig') & get(H.MaskedSpace_inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.MaskedSpace_inc_testnum,'Value',0);
end
set(H.MaskedSpace_FN,'String', FN.data);

% update stimulus durations
XStimParams.curr_stimdur2 = str2num(get(H.MaskedSpace_probeDUR, 'String'));
XStimParams.curr_stimdur = str2num(get(H.MaskedSpace_maskDUR, 'String'));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim