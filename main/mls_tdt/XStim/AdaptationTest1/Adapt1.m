% Adapt1: the GUI for Space tests of adaptation
% calc's a new stim on each rep and saves ACpower sequence to DATA file

%*******************************************************************************
%	The Graphical User Interface
%*******************************************************************************

XStimParams.reset_flag = 0;
GUI.locations1 = [];
GUI.locations2 = [];
% to start things off:
XStimParams.itdalone_flag = 0;
XStimParams.ildalone_flag = 0;
XStimParams.ABLalone_flag = 0;
XStimParams.space_flag = 1;
XStimParams.focalSPL = [0 30];

F.tot_ht = 0;
F.ind_ht = 20; 
%Figure window
if(~exist1('H.searchfig') &...
        ~exist1('H.ablfig') &...
        ~exist1('H.itdfig') &...
        ~exist1('H.ildfig') &...
        ~exist1('H.ildfreqfig') &...
        ~exist1('H.Adapt1fig'))
    temp = get(H.fig,'Position');
    F.space_l = temp(3) - 300; 
    F.space_wd = 400; 
    F.space_ht = 900;
    H.Adapt1fig = figure('Units','pixels',...
        'Position',[F.space_l 100 F.space_wd F.space_ht],...
        'Name','HRTF-based Space Measurement',...
        'NumberTitle','off');
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'ForegroundColor','blue',...
        'FontWeight','bold',...
        'String', 'Virtual Auditory Space Measurement');
    
    %reminder
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'ForegroundColor',[0 1 1],...
        'FontSize',12,...
        'FontWeight','bold',...
        'String', 'Spatial Filters in DSPs: use *.eq HRIRs');
    
    %Specify HRTF and Location files
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Specify Files');
    
    %Fully-cued HRTF file
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.spacefilepb = uicontrol('Parent',H.Adapt1fig,...
        'Style','pushbutton',...
        'Units','pixels',...
        'Position',[2 F.space_ht-F.tot_ht 195 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Select Fully-cued Filter File',...
        'Callback','setInfo_Adapt1');
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.spacefilehdr = uicontrol('Parent',H.Adapt1fig,...
        'Style','text',...
        'Position',[0 F.space_ht-F.tot_ht 120 F.ind_ht],...
        'String', 'Fully-cued HRTF File:');
    H.spacefile = uicontrol('Parent',H.Adapt1fig,...
        'Style','text',...
        'Position',[120 F.space_ht-F.tot_ht F.space_wd-120 F.ind_ht],...
        'HorizontalAlignment','left',...
        'ForegroundColor','blue',...
        'String', FN.space_std);
    
    F.tot_ht = F.tot_ht + F.ind_ht*2;
    %Set Stimulus Parameters
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Set Stimulus Parameters');
    F.tot_ht = F.tot_ht + F.ind_ht;    
    
    
    uicontrol('Parent',H.Adapt1fig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[2 F.space_ht-F.tot_ht 110 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Choose Stim Type');
    H.Adapt1_stim_type = uicontrol('Parent',H.Adapt1fig,...
        'Style','popup',...
        'Units','pixels',...
        'Position',[112 F.space_ht-F.tot_ht 110 F.ind_ht],...
        'BackgroundColor','white',...
        'String',...
        'Tone|Gammatone| | | |NarrowBand Noise| |BroadBand Noise| |BBN-LNN| ',...
        'Value',get(H.stim_type,'Value'),...
        'Callback','setInfo_Adapt1');
    
    H.Adapt1_stimFNpb = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[222 F.space_ht-F.tot_ht 20 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','black',...
        'FontWeight','bold',...
        'String', 'FN',...
        'Callback','setInfo_Adapt1');
    
    H.Adapt1_stimFN = uicontrol('Parent',H.Adapt1fig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[242 F.space_ht - F.tot_ht F.space_wd-215 F.ind_ht],...
        'ForegroundColor','blue',...
        'String', FN.stim);
    
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.Adapt1_DynamicRangepb = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[222 F.space_ht-F.tot_ht F.space_wd-222 F.ind_ht],...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'String', 'push for High Dynamic Range',...
        'Callback','setinfo_Adapt1');
    
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.Adapt1_DynamicRangeTxt = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[222 F.space_ht-F.tot_ht F.space_wd-222 F.ind_ht],...
        'BackgroundColor',[.8 .8 .8],...
        'String', 'SPL now 20dB below normal');
    
    % curr_freq
    H.Adapt1_freqtxt = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 125 F.ind_ht],...
        'String', 'center frequency (Hz):');
    H.Adapt1_freq = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[125 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.curr_freq),...
        'Callback','setInfo_Adapt1');
    
    H.Adapt1_bandlimitstxt = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 125 F.ind_ht],...
        'String', 'Bandlimits for NBNoise');
    H.Adapt1_bandlimits1 = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[125 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.bandlimits(1)),...
        'Callback','SetBandLimits');
    H.Adapt1_bandlimits2 = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[190 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.bandlimits(2)),...
        'Callback','SetBandLimits');
    
    %ABL (atten)
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
        'String', 'Atten (dB):');
    H.Adapt1_ABL = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.curr_ABL),...
        'Callback','setInfo_Adapt1');
    %DUR
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
        'String', 'Stim Dur (ms):');
    H.Adapt1_DUR = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.curr_stimdur),...
        'Callback','setInfo_Adapt1');
    %ISI
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
        'String', 'ISI (ms):');
    H.Adapt1_ISI = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.test_ISI));
    %NumReps
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
        'String', '# Reps:');
    H.Adapt1_numreps = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.numreps),...
        'Callback','setInfo_Adapt1');
    
    % modulate  sound1
    F.tot_ht = F.tot_ht + F.ind_ht*2;
    switch XStimParams.mod_type
        case 'Tone'
            mod_num(1) = 1;
        case 'Sq wave'
            mod_num(1) = 2;
        case 'LP Noise'
            mod_num(1) = 3;
        case 'File'
            mod_num(1) = 4;
        case 'None'
            mod_num(1) = 5;
        otherwise
            mod_num(1) = 5;
    end
    uicontrol(        'Style','text',...
        'Units','pixels',...
        'FontWeight','bold',...
        'Position',[0 F.space_ht - F.tot_ht 40 F.ind_ht],...
        'String', 'Mod ?');
    H.Adapt1_mod_type(1) = uicontrol(        'Style','popup',...
        'Units','pixels',...
        'Position',[40 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'BackgroundColor','White',...
        'String',...
        'Tone|Sq wave|LP Noise|File|None',...
        'Value',mod_num(1),...
        'Callback','setInfo_Adapt1');
    H.Adapt1_mod_txt(1) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'HorizontalAlignment','right',...
        'Position',[100 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'Freq:');
    H.Adapt1_mod_freq(1) = uicontrol(        'Style','edit',...
        'Units','pixels',...
        'Position',[130 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', num2str(XStimParams.mod_freq(1)),...
        'Callback','setInfo_Adapt1');
    
    H.Adapt1_mod_txtA(1) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'Position',[160 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'depth');
    H.Adapt1_mod_depth(1) = uicontrol(        'Style','edit',...
        'Units','pixels',...
        'Position',[190 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', num2str(XStimParams.mod_depth(1)),...
        'Callback','setInfo_Adapt1');
    
        % modulate  sound2
    switch XStimParams.mod_type2
        case 'Tone'
            mod_num(2) = 1;
        case 'Sq wave'
            mod_num(2) = 2;
        case 'LP Noise'
            mod_num(2) = 3;
        case 'File'
            mod_num(2) = 4;
        case 'None'
            mod_num(2) = 5;
        otherwise
            mod_num(2) = 5;
    end
    H.Adapt1_mod_type(2) = uicontrol(        'Style','popup',...
        'Units','pixels',...
        'Position',[220 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'BackgroundColor','White',...
        'String',...
        'Tone|Sq wave|LP Noise|File|None',...
        'Value',mod_num(2),...
        'Callback','setInfo_Adapt1');
    H.Adapt1_mod_txt(2) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'HorizontalAlignment','right',...
        'Position',[280 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'Freq:');
    H.Adapt1_mod_freq(2) = uicontrol(        'Style','edit',...
        'Units','pixels',...
        'Position',[310 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', num2str(XStimParams.mod_freq(2)),...
        'Callback','setInfo_Adapt1');
    
    H.Adapt1_mod_txtA(2) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'Position',[340 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'depth');
    H.Adapt1_mod_depth(2) = uicontrol(        'Style','edit',...
        'Units','pixels',...
        'Position',[370 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', num2str(XStimParams.mod_depth(2)),...
        'Callback','setInfo_Adapt1');    
    
    
    F.tot_ht = F.tot_ht + F.ind_ht;
 % stim 1   
    H.Adapt1_mod_txtB(1) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'Position',[160 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'phase');
    H.Adapt1_mod_phase(1) = uicontrol(        'Style','edit',...
        'Units','pixels',...
        'Position',[190 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', num2str(XStimParams.mod_phase(1)),...
        'Callback','setInfo_Adapt1');
    
    H.Adapt1_mod_pb(1) = uicontrol(        'Style','pushbutton',...
        'Units','pixels',...
        'Position',[100 F.space_ht-F.tot_ht 25 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','blue',...
        'String', 'FN',...
        'Visible','off',...
        'Callback','setInfo_Adapt1');
  % stim 2
        H.Adapt1_mod_txtB(2) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'Position',[340 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'phase');
    H.Adapt1_mod_phase(2) = uicontrol(        'Style','edit',...
        'Units','pixels',...
        'Position',[370 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', num2str(XStimParams.mod_phase(2)),...
        'Callback','setInfo_Adapt1');
    
    H.Adapt1_mod_pb(2) = uicontrol(        'Style','pushbutton',...
        'Units','pixels',...
        'Position',[270 F.space_ht-F.tot_ht 25 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','blue',...
        'String', 'FN',...
        'Visible','off',...
        'Callback','setInfo_Adapt1');

    
    
    
    
    
    
    %%%% SETup stimuli  
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'SetUp stimuli');
    
    %SPL range for epoch#1
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'BackgroundColor',[.8 .8 .8],...
        'ForegroundColor','blue',...
        'FontSize',10,...
        'String', 'Set relative SPL-ranges (below) between 1-31 dB.');
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'BackgroundColor',[.8 .8 .8],...
        'ForegroundColor','blue',...
        'FontSize',10,...
        'String', 'Actual SPLs will be attained with Atten (above).');
     F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'BackgroundColor',[.8 .8 .8],...
        'ForegroundColor',[.5 0 1],...
        'FontSize',10,...
        'String', 'stim: 00000-SPL1(1)-SPL2(1)-SPL1(1)-SPL2(2)...-00000...');
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'BackgroundColor',[.8 .8 .8],...
        'ForegroundColor',[.5 0 1],...
        'FontSize',10,...
        'String', '...SPL1(2)-SPL2(1)-SPL1(2)-SPL2(2)...-00000');
   
    
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 100 F.ind_ht],...
        'FontSize',8,...
        'String', 'epoch#1 SPLs');
    H.Adapt1_SPL1 = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[100 F.space_ht-F.tot_ht 100 F.ind_ht],...
        'String', num2str(XStimParams.SPL1),...
        'Callback','setInfo_Adapt1');
    %  duration                                             and epoch#2
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[210 F.space_ht-F.tot_ht 90 F.ind_ht],...
        'FontSize',8,...
        'String', 'epoch#2 SPLs');
    H.Adapt1_SPL2 = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[300 F.space_ht-F.tot_ht 100 F.ind_ht],...
        'String', num2str(XStimParams.SPL2),...
        'Callback','setInfo_Adapt1');
    
    
    % duration for epoch#1
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[20 F.space_ht-F.tot_ht 80 F.ind_ht],...
        'FontSize',8,...
        'String', 'duration (msec)');
    H.Adapt1_epoch_duration(1) = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[100 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.epoch_duration(1)),...
        'Callback','setInfo_Adapt1');
    %  duration                                             and epoch#2
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[210 F.space_ht-F.tot_ht 90 F.ind_ht],...
        'FontSize',8,...
        'String', 'duration (msec)');
    H.Adapt1_epoch_duration(2) = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[300 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.epoch_duration(2)),...
        'Callback','setInfo_Adapt1');
    
    % randomize stimuli?
    F.tot_ht = F.tot_ht + F.ind_ht; 
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[230 F.space_ht-F.tot_ht 100 F.ind_ht],...
        'BackGroundColor',[.8 .8 .8],...
        'ForeGroundColor','blue',...
        'String', 'randomize stims?');
    H.Adapt1_randStimFlag = uicontrol('Style','checkbox',...
        'Units','pixels',...
        'BackGroundColor',[.8 .8 .8],...
        'Position',[330 F.space_ht-F.tot_ht 20 F.ind_ht],...
        'value',0);

    
    
    % duration for silence periods
    F.tot_ht = F.tot_ht + F.ind_ht*.8;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[20 F.space_ht-F.tot_ht 200 F.ind_ht],...
        'FontSize',10,...
        'BackGroundColor',[.8 .8 .8],...
        'ForeGroundColor','blue',...
        'String', 'silence durations (msec): ');
    
       % plot stimulus in separate window
    H.Adapt1_plotstim = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[240 F.space_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','green',...
        'ForegroundColor','blue',...
        'String', 'plot stim',...
        'Callback','Adapt1_plotstim');

    
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'FontSize',8,...
        'String', 'lead');
    H.Adapt1_silence_dur(1) = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[60 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.silence_lead),...
        'Callback','setInfo_Adapt1');
    
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[130 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'FontSize',8,...
        'String', 'mid');
    H.Adapt1_silence_dur(2) = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[190 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.silence_mid),...
        'Callback','setInfo_Adapt1');
    
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[260 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'FontSize',8,...
        'String', 'trail');
    H.Adapt1_silence_dur(3) = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[320 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.silence_trail),...
        'Callback','setInfo_Adapt1');

    
    %%%% RUN Test    
    F.tot_ht = F.tot_ht + F.ind_ht*1.5;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Run Test');
    
    F.tot_ht = F.tot_ht + F.ind_ht;
    
    % increment test number
    H.Adapt1_inc_testnum = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[250 F.space_ht-F.tot_ht 140 F.ind_ht],...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'FontWeight','bold',...
        'String', 'INCREMENT testnum',...
        'Callback','setInfo_Adapt1');
    
    
    %Engage TEST
    H.engageAdapt1 = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','green',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'ENGAGE',...
        'Callback','Engage_Adapt1');
    
    %Record Data?
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[85 F.space_ht-F.tot_ht 100 F.ind_ht],...
        'String', 'Record Data?');
    H.Adapt1_recorddata = uicontrol('Style','checkbox',...
        'Units','pixels',...
        'Position',[185 F.space_ht-F.tot_ht 20 F.ind_ht],...
        'Callback','setInfo_Adapt1');
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.Adapt1_data_FN = uicontrol('Style','text',...
        'Units','pixels',...
        'ForegroundColor','blue',...
        'Position',[85 F.space_ht-F.tot_ht 100 F.ind_ht],...
        'String', FN.data);
    
    %Reset
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.resetAdapt1 = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[F.space_wd-100 F.space_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','cyan',...
        'ForegroundColor','black',...
        'FontWeight','bold',...
        'String', 'Reset',...
        'Callback','Reset_Adapt1');
    
    %Pause TEST
    F.ind_ht = 20; 
    F.tot_ht = F.tot_ht + F.ind_ht *2;
    H.pauseAdapt1 = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht*2],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Pause');
    
    %Exit
    H.exitAdapt1 = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[F.space_wd-100 F.space_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','yellow',...
        'FontWeight','bold',...
        'String', 'EXIT',...
        'Callback','Exit_Adapt1');
    
    %Status bar
    F.tot_ht = F.tot_ht + F.ind_ht *2;
    H.Adapt1_status = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Status: setting parameters');
    
    %Remaining Reps and Trials
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.space_ht-F.tot_ht 100 F.ind_ht],...
        'String', 'Remaining Reps:');
    H.Adapt1_remreps = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[102 F.space_ht-F.tot_ht 40 F.ind_ht],...
        'String', get(H.Adapt1_numreps,'String'));
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[144 F.space_ht-F.tot_ht 60 F.ind_ht],...
        'String', 'Trials:');
    H.Adapt1_remtrials = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[206 F.space_ht-F.tot_ht 40 F.ind_ht],...
        'String', '');
end %end GUI specification

% Initialize Application and get AP2 and XBUS locks
if(S232('S2init', 0, 'INIT_PRIMARY', 20000) == 0)
    disp('Cannot initialize a secondary process')
    return;
end

if(S232('APlock', 100, 0) == 0)
    disp('Cannot acquire lock on AP2 Card')
    s232('S2close');
    return;
end

if(S232('XBlock', 100, 0) == 0)
    disp('Cannot acquire lock on X Bus')
    s232('APunlock', 0);
    s232('S2close');
    return;
end

setInfo_Adapt1;