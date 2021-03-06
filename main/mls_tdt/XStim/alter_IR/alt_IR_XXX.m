%altIR: the GUI for altIR tests

%*******************************************************************************
%	The Graphical User Interface
%*******************************************************************************
% altered to use 12th octave bands

GUI.locations1 = [];
GUI.locations2 = [];
%cF = round(1000*exp(log(2)*([3:10]/3)));
cF = round(1000*exp(([12:40]/12)*log(2)))';
n_cF = length(cF);
XStimParams.reset_flag = 0;
% to start things off:
XStimParams.itdalone_flag = 0;
XStimParams.ildalone_flag = 0;
XStimParams.ABLalone_flag = 0;
XStimParams.space_flag = 1;

F.pos_incr = 20;
F.tot_ht = 0;
F.ind_ht = 20; 
%Figure window
if(~exist1('H.searchfig') &...
      ~exist1('H.ablfig') &...
      ~exist1('H.itdfig') &...
      ~exist1('H.ildfig') &...
      ~exist1('H.ildfreqfig') &...
      ~exist1('H.altIRfig'))
   temp = get(H.fig,'Position');
   F.altIR_l = temp(3) - 300; 
   F.altIR_wd = 400; 
   F.altIR_ht = 850;
   H.altIRfig = figure('Units','pixels',...
      'Position',[F.altIR_l 100 F.altIR_wd F.altIR_ht],...
      'Name','Space or ILDAlone Measurement with Alterred IRs',...
      'NumberTitle','off');

F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht F.altIR_wd F.ind_ht],...
   'ForegroundColor','blue',...
   'FontWeight','bold',...
   'String', 'Virtual Auditory Space Measurement with Alterred IRs');

%reminder
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht F.altIR_wd F.ind_ht],...
   'ForegroundColor','blue',...
   'String', 'Earphone Inverse Filters loaded to DSPs. Use *.STD HRTFs');

%Specify HRTF and Location files
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht F.altIR_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Specify Files');

%Fully-cued HRTF file
F.tot_ht = F.tot_ht + F.ind_ht;
H.altIR_spacefilehdr = uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Position',[0 F.altIR_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'Fully-cued HRTF File:');
H.altIR_spacefile = uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Position',[120 F.altIR_ht-F.tot_ht F.altIR_wd-120 F.ind_ht],...
   'HorizontalAlignment','left',...
   'ForegroundColor','blue',...
   'String', FN.space_std);
F.tot_ht = F.tot_ht + F.ind_ht;
H.altIR_spacefilepb = uicontrol('Parent',H.altIRfig,...
   'Style','pushbutton',...
   'Units','pixels',...
   'Position',[2 F.altIR_ht-F.tot_ht 195 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select Fully-cued Filter File',...
   'Callback','setInfo_altIR');

uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[200 F.altIR_ht-F.tot_ht 110 F.ind_ht],...
   'String', 'Fully-cued Test');
H.altIR_spaceflag = uicontrol('Parent',H.altIRfig,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[315 F.altIR_ht-F.tot_ht 30 F.ind_ht],...
   'Value',XStimParams.space_flag,...
   'Callback','setInfo_altIR');

% allow saving of altered HRTFs
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[2 F.altIR_ht-F.tot_ht 130 F.ind_ht],...
   'String', 'save altered HRTFs?');
H.altIR_saveHRTFflag = uicontrol('Parent',H.altIRfig,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[135 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
   'Value',XStimParams.space_flag,...
   'Callback','setInfo_altIR');
uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[200 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'Filename:');
H.altIR_saveHRTF_FN = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[280 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','white',...
   'String', FN.saveHRTF,...
   'Callback','setInfo_altIR');
% use saved altered HRTFs?
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[2 F.altIR_ht-F.tot_ht 130 F.ind_ht],...
   'String', 'use saved altered HRTFs?');
H.altIR_use_saveHRTFflag = uicontrol('Parent',H.altIRfig,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[135 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
   'Value',XStimParams.use_saveHRTFflag,...
   'Callback','setInfo_altIR');
% use saved alterations?
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[2 F.altIR_ht-F.tot_ht 130 F.ind_ht],...
   'String', 'use saved ALTerations?');
H.altIR_use_saveALTflag = uicontrol('Parent',H.altIRfig,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[135 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
   'Value',XStimParams.use_saveALTflag,...
   'Callback','setInfo_altIR');


%ILDAlone filter file
F.tot_ht = F.tot_ht + F.ind_ht + 10;
H.altIR_ildalonefilehdr = uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Position',[0 F.altIR_ht-F.tot_ht 120 F.ind_ht],...
    'ForegroundColor','blue',...
   'String', 'ILDAlone File:');
H.altIR_ildalonefile = uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Position',[120 F.altIR_ht-F.tot_ht F.altIR_wd-120 F.ind_ht],...
   'HorizontalAlignment','left',...
   'ForegroundColor','blue',...
   'String', FN.ildalone_std);
F.tot_ht = F.tot_ht + F.ind_ht;
H.altIR_ildalonepb = uicontrol('Parent',H.altIRfig,...
   'Style','pushbutton',...
   'Units','pixels',...
   'Position',[2 F.altIR_ht-F.tot_ht 195 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select ILDAlone Filter File',...
   'Callback','setinfo_altIR');

uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[200 F.altIR_ht-F.tot_ht 110 F.ind_ht],...
   'String', 'ILDAlone Test');
H.altIR_ildaloneflag = uicontrol('Parent',H.altIRfig,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[315 F.altIR_ht-F.tot_ht 30 F.ind_ht],...
   'Value',~(XStimParams.space_flag),...
   'Callback','setInfo_altIR');

% set stimulus params
F.tot_ht = F.tot_ht + F.ind_ht*2;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht F.altIR_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Set Stimulus Parameters');
F.tot_ht = F.tot_ht + F.ind_ht;

H.altIR_stimFNpb = uicontrol('Parent',H.altIRfig,...
   	'Style','pushbutton',...
   	'Units','pixels',...
   	'Position',[2 F.altIR_ht-F.tot_ht 195 F.ind_ht],...
   	'BackgroundColor','red',...
   	'ForegroundColor','white',...
   	'FontWeight','bold',...
   	'String', 'Select Stimulus File',...
   	'Callback','setinfo_altIR');
H.altIR_stimFN = uicontrol('Parent',H.altIRfig,...
   	'Style','text',...
   	'Units','pixels',...
   	'Position',[F.altIR_wd-175 F.altIR_ht - F.tot_ht 175 F.ind_ht],...
   	'ForegroundColor','blue',...
   	'String', FN.stim);

F.tot_ht = F.tot_ht + F.ind_ht*2;

H.altIR_NBstackpb = uicontrol('Style','pushbutton',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 150 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Bands for NBNoise Stack',...
   'Callback','NBNoiseStackBands');

F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 200 F.ind_ht],...
   'String', 'Bandwidth for BPNoise (octaves)');
H.altIR_bandwidth = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[205 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.bandwidth),...
   'Callback','SetBandWidth');
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 150 F.ind_ht],...
   'String', 'Bandlimits for NBNoise');
H.altIR_bandlimits1 = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[155 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.bandlimits(1)),...
   'Callback','SetBandLimits');
H.altIR_bandlimits2 = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[240 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.bandlimits(2)),...
   'Callback','SetBandLimits');

%ITD
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'ITD (us):');

H.altIR_ITD = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_ITD));

uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[145 F.altIR_ht-F.tot_ht 150 F.ind_ht],...
   'String', 'Hold ITD constant?');
H.altIR_ITDchk = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[300 F.altIR_ht-F.tot_ht 30 F.ind_ht],...
   'Value', XStimParams.space_flag);
%ABL
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'ABL (dB):');
H.altIR_ABL = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_ABL));
%DUR
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'Stim Dur (ms):');
H.altIR_DUR = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_stimdur));
%ISI
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'ISI (ms):');
H.altIR_ISI = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.test_ISI));
%NumReps
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'String', '# Reps:');
H.altIR_numreps = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.numreps),...
   'Callback','setinfo_altIR');


F.tot_ht = F.tot_ht + F.ind_ht*2;
% control IR alteration
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[200 F.altIR_ht-F.tot_ht 70 F.ind_ht],...
   'String', 'Alter IRs?');
uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[270 F.altIR_ht-F.tot_ht 40 F.ind_ht],...
   'String', 'left IR');
H.alter_left = uicontrol('Parent',H.altIRfig,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[310 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
   'Value',XStimParams.alter_left,...
   'Callback','setInfo_altIR');
uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[330 F.altIR_ht-F.tot_ht 40 F.ind_ht],...
   'String', 'right IR');
H.alter_right = uicontrol('Parent',H.altIRfig,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[370 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
   'Value',XStimParams.alter_right,...
   'Callback','setInfo_altIR');

if 0
uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[250 F.altIR_ht-F.tot_ht 100 F.ind_ht],...
   'String', 'set bandlimits?');
H.alter_setbands = uicontrol('Parent',H.altIRfig,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[350 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
   'Value',XStimParams.alter_setbands,...
   'Callback','setInfo_altIR');
end

% #12th octave bands to combine
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 100 F.ind_ht],...
   'String', '# 12ths to combine?');
H.altIR_Twelths_to_combine = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[100 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.Twelths_to_combine),...
   'Callback','setInfo_altIR');

% alter magnitude
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Parent',H.altIRfig,...
   'Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'String', 'magnitude');
H.alterMag = uicontrol('Parent',H.altIRfig,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[60 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
   'Value',XStimParams.alterMag,...
   'Callback','setInfo_altIR');
H.alterMag_dB = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[80 F.altIR_ht-F.tot_ht 40 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.alterMag_dB),...
   'Callback','setInfo_altIR');
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[120 F.altIR_ht-F.tot_ht 30 F.ind_ht],...
   'String', 'dB');

if 0
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[200 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'String', 'Bandlimits');
H.alterMag_Minbandlimits(1) = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[260 F.altIR_ht-F.tot_ht 70 F.ind_ht],...
   'BackgroundColor','white',...
   'String',XStimParams.alterMag_Minbandlimits(1),...
   'Callback','setInfo_altIR');
H.alterMag_Maxbandlimits(1) = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[330 F.altIR_ht-F.tot_ht 70 F.ind_ht],...
   'BackgroundColor','white',...
   'String',XStimParams.alterMag_Maxbandlimits(1),...
   'Callback','setInfo_altIR');
end
F.tot_ht = F.tot_ht + F.ind_ht;


for i = 1:n_cF-1
    H.altIR_Mag_octlist(i) = uicontrol('Parent',H.altIRfig,...
    'Style','checkbox',...
    'Units','pixels',...
    'Position',[(i-1)*14+2 F.altIR_ht-F.tot_ht 12 F.ind_ht],...
    'Value',XStimParams.Mag_octlistBIN(i),...
    'Callback','setInfo_altIR');
end
F.tot_ht = F.tot_ht + F.ind_ht;

for i = 1:4:n_cF-1
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[(i-1)*14 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
        'String', num2str(roundn(cF(i),50)/1000));
end

F.tot_ht = F.tot_ht + F.ind_ht;

% alter phase
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'String', 'phase');
H.alterPhase = uicontrol('Parent',H.altIRfig,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[60 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
   'Value',XStimParams.alterPhase,...
   'Callback','setInfo_altIR');

H.alterPhase_rad= uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[80 F.altIR_ht-F.tot_ht 50 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.alterPhase_rad),...
   'Callback','setInfo_altIR');
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[130 F.altIR_ht-F.tot_ht 40 F.ind_ht],...
   'String', 'radians');
if 0
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[200 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'String', 'Bandlimits');
H.alterPhase_Minbandlimits(1) = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[260 F.altIR_ht-F.tot_ht 70 F.ind_ht],...
   'BackgroundColor','white',...
   'String',XStimParams.alterMag_Minbandlimits(1),...
   'Callback','setInfo_altIR');
H.alterPhase_Maxbandlimits(1) = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[330 F.altIR_ht-F.tot_ht 70 F.ind_ht],...
   'BackgroundColor','white',...
   'String',XStimParams.alterMag_Maxbandlimits(1),...
   'Callback','setInfo_altIR');
end

F.tot_ht = F.tot_ht + F.ind_ht;
for i = 1:n_cF-1
    H.altIR_Phase_octlist(i) = uicontrol('Parent',H.altIRfig,...
    'Style','checkbox',...
    'Units','pixels',...
    'Position',[(i-1)*14+2 F.altIR_ht-F.tot_ht 12 F.ind_ht],...
    'Value',XStimParams.Phase_octlistBIN(i),...
    'Callback','setInfo_altIR');
end
F.tot_ht = F.tot_ht + F.ind_ht;

for i = 1:4:n_cF-1
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[(i-1)*14 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
        'String', num2str(roundn(cF(i),50)/1000));
end

% Only build files?? (do not play them out)
F.tot_ht = F.tot_ht + F.ind_ht*2;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[15 F.altIR_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'ONLY Build Files?');
H.altIR_buildOnly = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[135 F.altIR_ht-F.tot_ht 20 F.ind_ht],...
   'Callback','setInfo_altIR');

% run test
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht F.altIR_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Run Test');


%Engage TEST
F.tot_ht = F.tot_ht + F.ind_ht;
H.engagealtIR = uicontrol('Style','pushbutton',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','green',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'ENGAGE',...
   'Callback','Engage_altIR');

%Record Data?
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[85 F.altIR_ht-F.tot_ht 100 F.ind_ht],...
   'String', 'Record Data?');
H.altIR_recorddata = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[190 F.altIR_ht-F.tot_ht 30 F.ind_ht],...
   'Callback','setInfo_altIR');
% increment test number
H.altIR_inc_testnum = uicontrol('Style','pushbutton',...
   'Units','pixels',...
   'Position',[250 F.altIR_ht-F.tot_ht 140 F.ind_ht],...
   'BackgroundColor','yellow',...
   'ForegroundColor','blue',...
   'FontWeight','bold',...
   'String', 'INCREMENT testnum',...
   'Callback','setInfo_altIR');
%Record Data FN
F.tot_ht = F.tot_ht + F.ind_ht;
H.altIR_recorddata_FN = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[85 F.altIR_ht-F.tot_ht 100 F.ind_ht],...
   'ForegroundColor','blue',...
   'String', FN.data);


%Pause TEST
F.tot_ht = F.tot_ht + F.ind_ht;
H.pausealtIR = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Pause');

%Reset
H.resetaltIR = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[F.altIR_wd-100 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','cyan',...
   'ForegroundColor','black',...
   'FontWeight','bold',...
   'String', 'Reset',...
   'Callback','Reset_altIR');

%Build/Play Trials
F.tot_ht = F.tot_ht + F.ind_ht + F.pos_incr;
H.altIR_buildplay = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 175 F.ind_ht],...
   'BackgroundColor','yellow',...
   'ForegroundColor','blue',...
   'String', 'Build/Play Status');

%Exit
H.exitaltIR = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[F.altIR_wd-100 F.altIR_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','yellow',...
   'FontWeight','bold',...
   'String', 'EXIT',...
   'Callback','Exit_altIR');

%Remaining Reps and Trials
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 140 F.ind_ht],...
   'String', '    Remaining   reps:');
H.altIR_remreps = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[142 F.altIR_ht-F.tot_ht 40 F.ind_ht],...
   'String', get(H.altIR_numreps,'String'));
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[184 F.altIR_ht-F.tot_ht 60 F.ind_ht],...
   'String', '   Trials:');
H.altIR_remtrials = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[246 F.altIR_ht-F.tot_ht 40 F.ind_ht],...
   'String', '');

%Plotting
F.tot_ht = F.tot_ht + F.ind_ht;
H.altIR_plotpsd = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[0 F.altIR_ht-F.tot_ht 175 F.ind_ht],...
   'BackgroundColor','yellow',...
   'ForegroundColor','black',...
   'FontWeight','bold',...
   'String', 'Plot PSD');

end %end GUI specification

stimuli_dir = FN.temp_stim_path;

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

SetInfo_altIR;