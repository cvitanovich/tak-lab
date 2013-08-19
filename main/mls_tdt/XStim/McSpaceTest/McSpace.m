% McSpace: the GUI for Space tests using stimuli simuliar to Dean and McAlpine's
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
if(~exist1('H.searchfig') &&...
      ~exist1('H.ablfig') &&...
      ~exist1('H.itdfig') &&...
      ~exist1('H.ildfig') &&...
      ~exist1('H.ildfreqfig') &&...
      ~exist1('H.mcSpacefig'))
   temp = get(H.fig,'Position');
   F.space_l = temp(3) - 300; 
   F.space_wd = 400; 
   F.space_ht = 900;
   H.mcSpacefig = figure('Units','pixels',...
      'Position',[F.space_l 100 F.space_wd F.space_ht],...
      'Name','HRTF-based Space Measurement with McAlpine stimuli',...
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
H.mcSpace_filepb = uicontrol('Parent',H.mcSpacefig,...
   'Style','pushbutton',...
   'Units','pixels',...
   'Position',[2 F.space_ht-F.tot_ht 195 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select Fully-cued Filter File',...
   'Callback','setInfo_mcSpace');
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_filehdr = uicontrol('Parent',H.mcSpacefig,...
   'Style','text',...
   'Position',[0 F.space_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'Fully-cued HRTF File:');
H.mcSpace_file = uicontrol('Parent',H.mcSpacefig,...
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

uicontrol('Parent',H.mcSpacefig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[2 F.space_ht-F.tot_ht 110 F.ind_ht],...
   	'BackgroundColor','red',...
   	'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Choose Stim Type');
H.mcSpace_stim_type = uicontrol('Parent',H.mcSpacefig,...
   'Style','popup',...
   'Units','pixels',...
   'Position',[112 F.space_ht-F.tot_ht 110 F.ind_ht],...
   'BackgroundColor','white',...
   'String',...
   'Tone|Gammatone|Tone Stack|GammaTone Stack|BandPassed Noise|NarrowBand Noise|NarrowBand Noise Stack|BroadBand Noise|File|BBN - LNN|None',...
   'Value',get(H.stim_type,'Value'),...
   'Callback','setInfo_mcSpace');

H.mcSpace_stimFNpb = uicontrol('Style','pushbutton',...
   'Units','pixels',...
   'Position',[222 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','black',...
   'FontWeight','bold',...
   'String', 'FN',...
   'Callback','setInfo_mcSpace');

H.mcSpace_stimFN = uicontrol('Parent',H.mcSpacefig,...
   	'Style','text',...
   	'Units','pixels',...
   	'Position',[242 F.space_ht - F.tot_ht F.space_wd-215 F.ind_ht],...
    'ForegroundColor','blue',...
   	'String', FN.stim);

F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_repeatmodeTxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[5 F.space_ht-F.tot_ht 100 F.ind_ht],...
   	'BackgroundColor',[.8 .8 .8],...
   'String', 'same stim each rep?');
H.mcSpace_repeatmode = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[107 F.space_ht-F.tot_ht 15 F.ind_ht],...
   	'BackgroundColor',[.8 .8 .8],...
    'Value', 0, ...
   'Callback','setInfo_mcSpace');

H.mcSpace_DynamicRangepb = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[222 F.space_ht-F.tot_ht F.space_wd-222 F.ind_ht],...
   'BackgroundColor','yellow',...
   'ForegroundColor','blue',...
   'String', 'push for High Dynamic Range',...
   'Callback','setinfo_mcSpace');
F.tot_ht = F.tot_ht + F.ind_ht;

H.mcSpace_DynamicRangeTxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[222 F.space_ht-F.tot_ht F.space_wd-222 F.ind_ht],...
   	'BackgroundColor',[.8 .8 .8],...
   'String', 'SPL now 20dB below normal');

% curr_freq
H.mcSpace_freqtxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 125 F.ind_ht],...
   'String', 'center frequency (Hz):');
H.mcSpace_freq = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[125 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.curr_freq),...
   'Callback','setInfo_mcSpace');

F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_bwtxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 205 F.ind_ht],...
   'String', 'Bandwidth for BPNoise (octaves)');
H.mcSpace_bandwidth = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[205 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.bandwidth),...
   'Callback','setInfo_mcSpace');
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_bltxt(1) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 155 F.ind_ht],...
   'String', 'Bandlimits for NBNoise');
H.mcSpace_bandlimits(1) = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[155 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.bandlimits(1)),...
   'Callback','setInfo_mcSpace');
H.mcSpace_bandlimits(2) = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[240 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.bandlimits(2)),...
   'Callback','setInfo_mcSpace');

F.tot_ht = F.tot_ht + 3*F.ind_ht;
%ABL (atten)
H.mcSpace_ABLtxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'Atten (dB):');
H.mcSpace_ABL = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_ABL),...
   'Callback','setInfo_mcSpace');

%DUR
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_DURtxt(1) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'Stim Dur (ms):');
H.mcSpace_DUR = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_stimdur),...
   'Callback','setInfo_mcSpace');
H.mcSpace_DURtxt(2) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[142 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'BackgroundColor',[.8 .8 .8], ...
   'ForegroundColor',[.5 .5 1], ...
   'FontSize',10,...
   'String', '32K max');

% ramppts supplementary text
H.mcSpace_rampptstxt(2) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[200 F.space_ht-F.tot_ht 190 80],...
   'BackgroundColor',[.8 .8 .8], ...
   'ForegroundColor',[1 .5 .5], ...
   'FontSize',12,...
   'String', 'suppl text');

%ISI
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_ISItxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'ISI (ms):');
H.mcSpace_ISI = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.test_ISI));

% text for choosing McTestType
H.mcSpace_TestTypetxt(1) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[145 F.space_ht-F.tot_ht 255 F.ind_ht],...
   'String', '             AMdpth:       AMper:       ILD    Mono');
% text for choosing McTestType
H.mcSpace_TestTypetxt(2) = uicontrol('Style','text',...
   'Units','pixels',...
   'Fontsize',8,...
   'Position',[145 F.space_ht-F.tot_ht-F.ind_ht/2 255 F.ind_ht/1.5],...
   'String', 'SPL           lin  dB         lin dB                L  R ');

%NumReps
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_numrepstxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', '# Reps:');
H.mcSpace_numreps = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.numreps),...
   'Callback','setInfo_mcSpace');

F.tot_ht = F.tot_ht + F.ind_ht*.5;
H.mcSpace_TestTypetxt(3) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[145 F.space_ht-F.tot_ht 25 F.ind_ht],...
   'String', 'Test:');
 
H.mcSpace_McTestType(1) = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[170 F.space_ht-F.tot_ht 15 F.ind_ht],...
   'Callback','setInfo_mcSpace');
H.mcSpace_McTestType(2) = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[220 F.space_ht-F.tot_ht 15 F.ind_ht],...
   'Callback','setInfo_mcSpace');
H.mcSpace_McTestType(3) = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[235 F.space_ht-F.tot_ht 15 F.ind_ht],...
   'Callback','setInfo_mcSpace');
H.mcSpace_McTestType(4) = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[275 F.space_ht-F.tot_ht 15 F.ind_ht],...
   'Callback','setInfo_mcSpace');
H.mcSpace_McTestType(5) = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[290 F.space_ht-F.tot_ht 15 F.ind_ht],...
   'Callback','setInfo_mcSpace');
H.mcSpace_McTestType(6) = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[325 F.space_ht-F.tot_ht 15 F.ind_ht],...
   'Callback','setInfo_mcSpace');
H.mcSpace_McTestType(7) = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[350 F.space_ht-F.tot_ht 15 F.ind_ht],...
   'Callback','setInfo_mcSpace');
H.mcSpace_McTestType(8) = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[365 F.space_ht-F.tot_ht 15 F.ind_ht],...
   'Callback','setInfo_mcSpace');

F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_TestTypetxt(4) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[220 F.space_ht-F.tot_ht 125 F.ind_ht],...
   'String', '#epochs b/wn probes');
H.mcSpace_offEpochs = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[345 F.space_ht-F.tot_ht 55 F.ind_ht],...
   'String', num2str(XStimParams.offEpochs),...
   'Callback','setInfo_mcSpace');

% modulate  sound
    F.tot_ht = F.tot_ht + F.ind_ht*1.5;
    switch XStimParams.mod_type
        case 'Tone'
            mod_num = 1;
        case 'Sq wave'
            mod_num = 2;
        case 'LP Noise'
            mod_num = 3;
        case 'File'
            mod_num = 4;
        case 'None'
            mod_num = 5;
        otherwise
            mod_num = 5;
    end
    uicontrol(        'Style','text',...
        'Units','pixels',...
        'FontWeight','bold',...
        'Position',[0 F.space_ht - F.tot_ht 80 F.ind_ht],...
        'String', '   Modulation:');
    H.mcSpace_mod_type = uicontrol(        'Style','popup',...
        'Units','pixels',...
        'Position',[80 F.space_ht-F.tot_ht 100 F.ind_ht],...
        'BackgroundColor','White',...
        'String',...
        'Tone|Sq wave|LP Noise|File|None',...
        'Value',mod_num,...
        'Callback','setInfo_mcSpace');
    H.mcSpace_mod_txt(1) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'HorizontalAlignment','right',...
        'Position',[180 F.space_ht-F.tot_ht 40 F.ind_ht],...
        'String', 'Freq:');
    H.mcSpace_mod_freq = uicontrol(        'Style','edit',...
        'Units','pixels',...
        'Position',[220 F.space_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_freq(1)),...
        'Callback','setInfo_mcSpace');
    
    H.mcSpace_mod_txt(2) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'Position',[260 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'depth');
    H.mcSpace_mod_depth = uicontrol(        'Style','edit',...
        'Units','pixels',...
        'Position',[290 F.space_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_depth(1)),...
        'Callback','setInfo_mcSpace');

    F.tot_ht = F.tot_ht + F.ind_ht;
    H.mcSpace_mod_txt(3) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'Position',[10 F.space_ht-F.tot_ht 200 F.ind_ht],...
        'BackgroundColor',[1 .5 .5], ...
        'ForegroundColor',[.4 .4 1], ...
        'FontSize',12,...
        'Visible','off',...
        'String', 'use BBN - LNN ??');

    H.mcSpace_mod_txt(4) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'Position',[260 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'phase');
    H.mcSpace_mod_phase = uicontrol(        'Style','edit',...
        'Units','pixels',...
        'Position',[290 F.space_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_phase(1)),...
        'Callback','setInfo_mcSpace');

    H.mcSpace_mod_pb = uicontrol(        'Style','pushbutton',...
        'Units','pixels',...
        'Position',[200 F.space_ht-F.tot_ht 25 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','blue',...
        'String', 'FN',...
        'Visible','off',...
        'Callback','setInfo_mcSpace');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  %%%% SETup SPLs  
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_SPLtxt(1) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'SetUp SPLs');
  
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_SPLtxt(2) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 140 F.ind_ht],...
   'FontSize',10,...
   'String', 'focal-range SPL');


H.mcSpace_SPLtxt(3) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[300 F.space_ht-F.tot_ht 90 F.ind_ht],...
   'FontSize',10,...
   'String', 'AM freq');

%SPL FOCAL-range lower-SPL slider   (slider#1)
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_SPLslider(1) = uicontrol('Style','slider',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd-60 F.ind_ht *.75],...
   'BackgroundColor','blue',...
   'Min',0,'Max',30,...
   'SliderStep',[1/31 5/31],...
   'Value',XStimParams.focalSPL(1),...
   'Callback','setInfo_mcSpace');
H.mcSpace_SPLslider_txt(1) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[F.space_wd-60 F.space_ht-F.tot_ht 60 F.ind_ht*.75],...
   'FontSize',8,...
   'String', [num2str(XStimParams.focalSPL(1)) ' dB']);
%SPL FOCAL-range upper-SPL slider  (slider#2)
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_SPLslider(2) = uicontrol('Style','slider',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd-60 F.ind_ht*.75],...
   'BackgroundColor','red',...
   'Min',0,'Max',30,...
   'SliderStep',[1/31 5/31],...
   'Value',XStimParams.focalSPL(2),...
   'Callback','setInfo_mcSpace');
H.mcSpace_SPLslider_txt(2) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[F.space_wd-60 F.space_ht-F.tot_ht 60 F.ind_ht*.75],...
   'FontSize',8,...
   'String', [num2str(XStimParams.focalSPL(2)) ' dB']);
    
  % FOCAL total probability
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_focalProbtxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 150 F.ind_ht],...
   'String', 'FOCAL total probability');
H.mcSpace_focalProb = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[155 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', num2str(XStimParams.focalProb),...
   'Callback','setInfo_mcSpace');

% epoch duration
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_epoch_durationtxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 150 F.ind_ht],...
   'String', 'epoch duration (msec)');
H.mcSpace_epoch_duration = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[155 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', num2str(XStimParams.epoch_duration(1)),...
   'Callback','setInfo_mcSpace');

% ramppts (added Jun 21, 2007)
H.mcSpace_rampptstxt(1) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[250 F.space_ht-F.tot_ht 70 F.ind_ht],...
   'String', 'ramppts : ');
H.mcSpace_ramppts = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[330 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.ramppts),...
   'Callback','setInfo_mcSpace');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
H.mcSpace_inc_testnum = uicontrol('Style','pushbutton',...
   'Units','pixels',...
   'Position',[250 F.space_ht-F.tot_ht 140 F.ind_ht],...
   'BackgroundColor','yellow',...
   'ForegroundColor','blue',...
   'FontWeight','bold',...
   'String', 'INCREMENT testnum',...
   'Callback','setInfo_mcSpace');

%Engage TEST
H.engage_mcSpace = uicontrol('Style','pushbutton',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','green',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'ENGAGE',...
   'Callback','Engage_mcSpace');

%Record Data?
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[85 F.space_ht-F.tot_ht 100 F.ind_ht],...
   'String', 'Record Data?');
H.mcSpace_recorddata = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[185 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'Callback','setInfo_mcSpace');
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_data_FN = uicontrol('Style','text',...
   'Units','pixels',...
   'ForegroundColor','blue',...
   'Position',[85 F.space_ht-F.tot_ht 100 F.ind_ht],...
   'String', FN.data);

%Reset
F.tot_ht = F.tot_ht + F.ind_ht;
H.resetmcSpace = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[F.space_wd-100 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','cyan',...
   'ForegroundColor','black',...
   'FontWeight','bold',...
   'String', 'Reset',...
   'Callback','Reset_mcSpace');

%Pause TEST
F.ind_ht = 20; 
F.tot_ht = F.tot_ht + F.ind_ht *2;
H.pausemcSpace = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht*2],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Pause');

%Exit
H.exitmcSpace = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[F.space_wd-100 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','yellow',...
   'FontWeight','bold',...
   'String', 'EXIT',...
   'Callback','Exit_mcSpace');

%Status bar
F.tot_ht = F.tot_ht + F.ind_ht *2;
H.mcSpace_status = uicontrol('Style','text',...
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
H.mcSpace_remreps = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[102 F.space_ht-F.tot_ht 40 F.ind_ht],...
   'String', get(H.mcSpace_numreps,'String'));
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[144 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', 'Trials:');
H.mcSpace_remtrials = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[206 F.space_ht-F.tot_ht 40 F.ind_ht],...
   'String', '');

% analyze cycle-by-cycle Vstrength
F.tot_ht = F.tot_ht + F.ind_ht;
H.mcSpace_analVstr = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[F.space_wd-100 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','yellow',...
   'FontWeight','bold',...
   'String', 'Analyze Vstr',...
   'Callback','anal_cycleVstr');

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

setInfo_mcSpace;