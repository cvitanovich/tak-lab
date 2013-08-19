%setInfo_LRsounds

XStimParams.stim_type = 'File';
XStimParams.test_type = 'LRsounds';

% engage stimbuilder
if get(H.LR_stimbuilderpb,'value')
    addpath('e:\kip\matlab\stimuli\stim_builder');
    addpath('e:\kip\matlab\SpacePicker2');
    stimbuilder;
    set(H.LR_stimbuilderpb,'value',0);
end


% pick stimulus FNs
if get(H.LR_stimFNpb(1),'value')
    [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.*'],'Select Left Ear stim File');
    set(H.LR_stimFN(1),'string', FN.stim);
end
D = dir([FN.stim_path FN.stim]);
dur(1) = round(1000 * (D.bytes/4) / TDT.Fs);
if get(H.LR_stimFNpb(2),'value')
    [FN.stim2,FN.stim_path2] = uigetfile([FN.stim_path2 '*.*'],'Select Right Ear stim File');
    set(H.LR_stimFN(2),'string', FN.stim2);
end
D = dir([FN.stim_path2 FN.stim2]);
dur(2) = round(1000 * (D.bytes/4) / TDT.Fs);
if any(diff(dur))
    disp(['sound durations do not all match - you may need to re-pick : ']);
    dur
    set(H.DUR,'foregroundcolor','red');
    XStimParams.curr_stimdur = [];
    set(H.DUR, 'String', num2str(dur));
else
    XStimParams.curr_stimdur = dur(1);
    set(H.DUR, 'String', num2str(XStimParams.curr_stimdur));
    set(H.DUR,'foregroundcolor','black');
end
clear D

%%%% get HRTF file
% check for fully-cued pushbutton (file-picker)
if get(H.LR_hrtfFNpb,'Value')
    addpath('e:\kip\matlab\SpacePicker2');
    if ~exist1('GUI2')    Globals_GUI2;   end
    [GUI2.hrtfFN,FN.hrtfPATH] = uigetfile([GUI2.hrtfPATH '*.*'],'Select HRTF File (*.std preferred)');
    if(GUI2.hrtfPATH ~= 0)
        set(H.LR_hrtfFN,'String',GUI2.hrtfFN);
    end
    set(H.LR_hrtfFNpb,'Value',0);

    if findstr(GUI2.hrtfFN,'std')
        set(H.LR_hrtftxt(1),'String', 'using ephone file')
    else
        set(H.LR_hrtftxt(1),'String', [])
    end
end

% engage spacepicker
if get(H.LR_spacepickerpb,'value')
    addpath('e:\kip\matlab\SpacePicker2');
    if ~exist1('GUI2')    Globals_GUI2;   end
    spacePicker2;
    set(H.LR_spacepickerpb,'value',0);
    set(H.LR_setITD,'value',1)
    set(H.LR_setILD,'value',1)
end

% ITD
XStimParams.curr_ITD = str2num(get(H.ITD,'string'));
low_itd  = max1([-250 str2num(get(H.loitd,'String'))]);
high_itd = min1([250 str2num(get(H.hiitd,'String'))]);
set(H.hiitd,'String',num2str(high_itd));
set(H.loitd,'String',num2str(low_itd));
XStimParams.hiitd = high_itd;
XStimParams.loitd = low_itd;
XStimParams.numitds = str2num(get(H.numitds,'String'));
step_itd = (XStimParams.hiitd - XStimParams.loitd)/(XStimParams.numitds-1);
if get(H.LR_setITD,'value')
    set(H.stepitd,'string',[]);
    XStimParams.itds = XStimParams.curr_ITD;
else    
    set(H.stepitd,'string',num2str(step_itd));
    XStimParams.itds = round(XStimParams.loitd:step_itd:XStimParams.hiitd);
end
clear low_itd high_itd step_itd

% ILD
XStimParams.curr_ILD = str2num(get(H.ILD,'string'));
low_ild  = max1([-20 str2num(get(H.loild,'String'))]);
high_ild = min1([20 str2num(get(H.hiild,'String'))]);
set(H.hiild,'String',num2str(high_ild));
set(H.loild,'String',num2str(low_ild));
XStimParams.hiild = high_ild;
XStimParams.loild = low_ild;
XStimParams.numilds = str2num(get(H.numilds,'String'));
step_ild = (XStimParams.hiild - XStimParams.loild)/(XStimParams.numilds-1);
if get(H.LR_setILD,'value')
    set(H.stepild,'string',[]);
    XStimParams.ilds = XStimParams.curr_ILD;
else    
    set(H.stepild,'string',num2str(step_ild));
    XStimParams.ilds = round(XStimParams.loild:step_ild:XStimParams.hiild);
end
clear low_ild high_ild step_ild

% ABL
ABL = max1([-120 min1([-20 str2num(get(H.ABL,'string'))])]);
XStimParams.curr_ABL = ABL;
set(H.ABL,'string',num2str(ABL))
low_abl  = max1([-120 str2num(get(H.loabl,'String'))]);
high_abl = min1([-20 str2num(get(H.hiabl,'String'))]);
set(H.hiabl,'String',num2str(high_abl));
set(H.loabl,'String',num2str(low_abl));
XStimParams.hiabl = high_abl;
XStimParams.loabl = low_abl;
XStimParams.numabls = str2num(get(H.numabls,'String'));   
step_abl = round((high_abl-low_abl)/(XStimParams.numabls-1));
if get(H.LR_setABL,'value')
    set(H.stepabl,'string',[]);
    XStimParams.abls = XStimParams.curr_ABL;
else    
    set(H.stepabl,'string',num2str(step_abl));
    XStimParams.abls = round(XStimParams.loabl:step_abl:XStimParams.hiabl);
end
clear low_abl high_abl step_abl

% ISI
XStimParams.test_ISI = str2num(get(H.ISI,'string'));

% numreps
XStimParams.numreps = str2num(get(H.numreps,'string'));

% increment test number
if(exist1('H.LRsoundsfig') & get(H.inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.inc_testnum,'Value',0);
end

% reset dataFN on this and xstim controlpanels
set(H.recorddata_FN, 'String', FN.data);
set(H.xstim_recorddata_FN, 'String', FN.data);

% save current parameter values
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])