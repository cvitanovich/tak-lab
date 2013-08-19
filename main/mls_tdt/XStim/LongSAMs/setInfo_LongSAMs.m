% SetInfo_LongSAMs
% combines all the other 'set...' scripts into one

% uses .eq HRIRs, therefore make sure earphone filters are not selected
XStimParams.ephone_flag = 0;

if 0
if get(H.LongSAMs_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.LongSAMs_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.LongSAMs_DynamicRangepb,'BackgroundColor','red');
    set(H.LongSAMs_DynamicRangeTxt,'visible','on');
else
    XStimParams.HiDynamicRange = 0;
    set(H.LongSAMs_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.LongSAMs_DynamicRangepb,'BackgroundColor','yellow');
    set(H.LongSAMs_DynamicRangeTxt,'visible','off');
end
end

XStimParams.buildOnly = 0;
XStimParams.space_flag = 1;
XStimParams.ildalone_flag = 0;
XStimParams.itdalone_flag = 0;
XStimParams.ABLalone_flag = 0;

% control visibility of file picking for ILA, ITA or fully-cued
set(H.LongSAMsfilehdr,'Enable','on');
set(H.LongSAMsfile,'Enable','on');
set(H.LongSAMsfilepb,'Enable','on');   

% check for fully-cued pushbutton (file-picker)
if isempty(FN.space_eq) | get(H.LongSAMsfilepb,'Value')
    [FN.space_eq,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File (*.eq preferred)');
end
if(FN.space_path ~= 0)
    set(H.LongSAMsfile,'String',FN.space_eq);
end
set(H.LongSAMsfilepb,'Value',0);
FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);

% set numreps
XStimParams.numreps = str2num(get(H.LongSAMs_numreps,'String'));

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
if(exist1('H.LongSAMsfig') & get(H.LongSAMs_inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.LongSAMs_inc_testnum,'Value',0);
end
set(H.LongSAMsdata_FN,'String', FN.data);

% update display of stim type
stim_val = 8;
XStimParams.stim_type = 'BroadBandNoise';

% update test_ISI
XStimParams.test_ISI = str2num(get(H.LongSAMs_ISI,'string'));

% check out modulation parameters for stim
if get(H.LongSAMstestmod_1or2, 'value') == 1        % doubly modulated
    XStimParams.LongSAMs_flag = 2;
    %set(H.LongSAMstestmod_1or2, 'string', 'push for single mod');
    set(H.LongSAMstestmod_txt(1),'visible','on')
    set(H.LongSAMstestmod_txt(4),'visible','on')
    set(H.LongSAMstestmod_freq(1),'visible','on')
    XStimParams.mod_freq(1) = str2num(get(H.LongSAMstestmod_freq(1),'String'));
    XStimParams.mod_freq(2) = str2num(get(H.LongSAMstestmod_freq(2),'string'));
elseif get(H.LongSAMstestmod_1or2, 'value') == 2    % singly modulated
    XStimParams.LongSAMs_flag = 1;
    XStimParams.mod_type = 'None';
    %set(H.LongSAMstestmod_1or2, 'string', 'push for double mod');
    set(H.LongSAMstestmod_txt(1),'visible','off')
    set(H.LongSAMstestmod_txt(4),'visible','off')
    set(H.LongSAMstestmod_freq(1),'visible','off')
    XStimParams.mod_freq(2) = str2num(get(H.LongSAMstestmod_freq(2),'string'));
elseif get(H.LongSAMstestmod_1or2, 'value') == 3    % preset quad-mod
    XStimParams.LongSAMs_flag = 3;
    set(H.LongSAMstestmod_txt(1),'visible','off')
    set(H.LongSAMstestmod_txt(4),'visible','off')
    set(H.LongSAMstestmod_freq(1),'visible','off')
    XStimParams.mod_freq(2) = str2num(get(H.LongSAMstestmod_freq(2),'string'));
    %XStimParams.mod_freq(2) = 1/16;
    %set(H.LongSAMstestmod_freq(2),'string',num2str(1/16))
    set(H.LongSAMstestmod_depth(1),'String',num2str(.75));
elseif get(H.LongSAMstestmod_1or2, 'value') == 4    % preset seven-mod
    XStimParams.LongSAMs_flag = 4;
    set(H.LongSAMstestmod_txt(1),'visible','off')
    set(H.LongSAMstestmod_txt(4),'visible','off')
    set(H.LongSAMstestmod_freq(1),'visible','off')
    XStimParams.mod_freq(2) = str2num(get(H.LongSAMstestmod_freq(2),'string'));
    %XStimParams.mod_freq(2) = 1/16;
    %set(H.LongSAMstestmod_freq(2),'string',num2str(1/16))
    set(H.LongSAMstestmod_depth(1),'String',num2str(.75));
end

% main (slow) modulation
XStimParams.mod_depth(1) = str2num(get(H.LongSAMstestmod_depth(1),'String'));
if ~get(H.LongSAMstestmod_linlog(1),'value')        % linear mod scale
    XStimParams.mod_type = 'Tone';
    set(H.LongSAMstestmod_linlog(1),'string','linear')
    set(H.LongSAMstestmod_txt(6),'visible','off')
    set(H.LongSAMstestmod_dB(1),'visible','off')
else                                            % dB mod scale
    XStimParams.mod_type = 'dBTone';
    set(H.LongSAMstestmod_linlog(1),'string','dB scale')
    set(H.LongSAMstestmod_txt(6),'visible','on')
    set(H.LongSAMstestmod_dB(1),'visible','on')
    XStimParams.mod_dB(1) = str2num(get(H.LongSAMstestmod_dB(1),'string'));
end    

%XStimParams.mod_freq(2) = str2num(get(H.LongSAMstestmod_freq(2),'String'));
tempStruct.buf_pts = round(TDT.Fs / XStimParams.mod_freq(2));
tempStruct.buf_pts = max1([min1([tempStruct.buf_pts 24 * TDT.Fs]) .01 * TDT.Fs]);

% test if freq#1 ok
if XStimParams.mod_freq(1) < XStimParams.mod_freq(2)
    XStimParams.mod_freq(1) = XStimParams.mod_freq(2);
end
XStimParams.mod_freq(1) = XStimParams.mod_freq(2) * round(XStimParams.mod_freq(1)/XStimParams.mod_freq(2));
temp = XStimParams.mod_freq(1) / XStimParams.mod_freq(2);

XStimParams.mod_freq(2) = TDT.Fs / tempStruct.buf_pts;
set(H.LongSAMstestmod_freq(2),'String',num2str(XStimParams.mod_freq(2)));

% test again if freq#1 ok
XStimParams.mod_freq(1) = XStimParams.mod_freq(2) * temp;
set(H.LongSAMstestmod_freq(1),'String',num2str(XStimParams.mod_freq(1)));



% update stimulus duration
tempStruct.nCycles_perRep = str2num(get(H.LongSAMstestmod_nCycles_perRep,'string'));
tempStruct.nptsTotalPlay = tempStruct.nCycles_perRep * tempStruct.buf_pts;
XStimParams.curr_stimdur = 1000 * tempStruct.nptsTotalPlay / TDT.Fs;
set(H.LongSAMs_DUR, 'String',num2str(XStimParams.curr_stimdur));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim