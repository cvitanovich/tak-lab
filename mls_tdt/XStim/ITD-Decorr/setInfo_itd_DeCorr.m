%setInfo_itd_decorr

XStimParams.stim_type = 'File';
XStimParams.test_type = 'ITD_decorr';
% increment test number
if(exist1('H.itd_decorrfig') & get(H.inc_testnum,'Value'))
	XStimParams.testnum = XStimParams.testnum +1;
	set(H.testnum, 'String', num2str(XStimParams.testnum))
   update_dataFN;
   set(H.inc_testnum,'Value',0);
end

% make list of possible test names (by columns)
TestNames{1} = 'silence';
TestNames{2} = 'Cohere loFrq';
TestNames{3} = 'Lvl1 loFrq';
TestNames{4} = 'Lvl2 LoFrq';
TestNames{5} = 'Cohere hiFrq';
TestNames{6} = 'Cohere hiFrq & Cohere LoFrq';
TestNames{7} = 'Cohere hiFrq & lvl1 LoFrq';
TestNames{8} = 'Cohere hiFrq & lvl2 LoFrq';
TestNames{9} = 'lvl1 hiFrq';
TestNames{10} = 'lvl1 hiFrq & Cohere LoFrq';
TestNames{11} = 'lvl1 hiFrq & lvl1 LoFrq';
TestNames{12} = 'lvl1 hiFrq & lvl2 LoFrq';
TestNames{13} = 'lvl2 hiFrq';
TestNames{14} = 'lvl2 hiFrq & Cohere LoFrq';
TestNames{15} = 'lvl2 hiFrq & lvl1 LoFrq';
TestNames{16} = 'lvl2 hiFrq & lvl2 LoFrq';

% check which stims to play
mrkr_clr{1}= 'blue';
mrkr_clr{2} = 'green';
mrkr_clr{3} = 'cyan';
mrkr_clr{4} = 'magenta';
mrkr_clr{5} = 'yellow';
mrkr_clr{6} = 'white';
mrkr_clr{7}= 'blue';
mrkr_clr{8} = 'green';
mrkr_clr{9} = 'cyan';
mrkr_clr{10} = 'magenta';
mrkr_clr{11} = 'yellow';
mrkr_clr{12} = 'white';

itest = 0;
XStimParams.legendstr = [];
for j = 1:4
    for i = 1:4
        XStimParams.stims_to_play(i,j) = get(H.stims_to_play(i,j),'value');
        if get(H.stims_to_play(i,j),'value')
            itest = itest+1;
            XStimParams.legendstr{itest} = TestNames{i+(j-1)*4};
            set(H.stims_to_play(i,j),'ForeGroundColor',mrkr_clr{itest})
            set(H.stims_to_play(i,j),'string', 'X')
        else
            set(H.stims_to_play(i,j),'string',[])
        end
    end
end
clear mrkr_clr TestNames i j itest

% set up noise FNs
dur = [XStimParams.curr_stimdur XStimParams.curr_stimdur];
if get(H.stim_FN4pb,'value')
    [FN.stim4L,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select coherent HiF stim File');
    FN.stim4R = FN.stim4L;
    D = dir([FN.stim_path FN.stim4L]);
    dur(1) = round(1000 * (D.bytes/4) / TDT.Fs);
    dur(2) = dur(1);
    set(H.stim_FN4,'string', FN.stim4L(1:end-4));
end
if get(H.stim_FN5pb,'value')
    [FN.stim5L,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select LEFT level1 HiF stim File');
    [FN.stim5R,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select RIGHT level1 HiF stim File');
    D = dir([FN.stim_path FN.stim5L]);
    dur(3) = round(1000 * (D.bytes/4) / TDT.Fs);
    D = dir([FN.stim_path FN.stim5R]);
    dur(4) = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.stim_FN5,'string', FN.stim5L(1:end-4));
end
if get(H.stim_FN6pb,'value')
    [FN.stim6L,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select LEFT level2 HiF stim File');
    [FN.stim6R,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select RIGHT level2 HiF stim File');
    D = dir([FN.stim_path FN.stim6L]);
    dur(5) = round(1000 * (D.bytes/4) / TDT.Fs);
    D = dir([FN.stim_path FN.stim6R]);
    dur(6) = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.stim_FN6,'string', FN.stim6L(1:end-4));
end
if get(H.stim_FN7pb,'value')
    [FN.stim7L,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select coherent LoF stim File');
    FN.stim7R = FN.stim7L;
    D = dir([FN.stim_path FN.stim7L]);
    dur(7) = round(1000 * (D.bytes/4) / TDT.Fs);
    D = dir([FN.stim_path FN.stim7R]);
    dur(8) = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.stim_FN7,'string', FN.stim7L(1:end-4));
end
if get(H.stim_FN8pb,'value')
    [FN.stim8L,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select LEFT level1 LoF stim File');
    [FN.stim8R,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select RIGHT level1 LoF stim File');
    D = dir([FN.stim_path FN.stim8L]);
    dur(9) = round(1000 * (D.bytes/4) / TDT.Fs);
    D = dir([FN.stim_path FN.stim8R]);
    dur(10) = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.stim_FN8,'string', FN.stim8L(1:end-4));
end
if get(H.stim_FN9pb,'value')
    [FN.stim9L,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select LEFT level2 LoF stim File');
    [FN.stim9R,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select RIGHT level2 LoF stim File');
    D = dir([FN.stim_path FN.stim9L]);
    dur(11) = round(1000 * (D.bytes/4) / TDT.Fs);
    D = dir([FN.stim_path FN.stim9R]);
    dur(12) = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.stim_FN9,'string', FN.stim9L(1:end-4));
end
if any(diff(dur))
    clc
    disp(['sound durations do not all match - you may need to re-pick : ']);
    disp(['#4) coherent HiFreq:     ' FN.stim4L '    '  FN.stim4R]);
    disp(['#5) level1 HiFreq:       ' FN.stim5L '    '  FN.stim5R]);
    disp(['#6) level2 HiFreq:       ' FN.stim6L '    '  FN.stim6R]);
    disp(['#7) coherent LoFreq:     ' FN.stim7L '    '  FN.stim7R]);
    disp(['#8) level1 LoFreq:       ' FN.stim8L '    '  FN.stim8R]);
    disp(['#9) level2 LoFreq:       ' FN.stim9L '    '  FN.stim9R]);
else
    if ~isempty(findstr(FN.stim_path,'partial_corr\'))
        D = dir([FN.stim_path FN.stim4L]);
    else
        D = dir([FN.stim_path 'partial_corr\' FN.stim4L]);
    end    
    XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.DUR, 'String', num2str(XStimParams.curr_stimdur));
end
clear D dur

%setITDstep
low_itd  = max1([-250 str2num(get(H.lowitd,'String'))]);
high_itd = min1([250 str2num(get(H.highitd,'String'))]);
set(H.highitd,'String',num2str(high_itd));
set(H.lowitd,'String',num2str(low_itd));
XStimParams.hiitd = high_itd;
XStimParams.loitd = low_itd;
XStimParams.numitds = str2num(get(H.numitds,'String'));
set(H.stepitd,'string',num2str((XStimParams.hiitd - XStimParams.loitd)/(XStimParams.numitds-1)));
clear low_itd high_itd

% reset dataFN on this and xstim controlpanels
set(H.recorddata_FN, 'String', FN.data);
set(H.xstim_recorddata_FN, 'String', FN.data);

% save current parameter values
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim