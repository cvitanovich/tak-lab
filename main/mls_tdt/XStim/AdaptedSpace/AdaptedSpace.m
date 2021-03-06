% AdaptedSpace: the GUI for AdaptedSpace tests
% For randomized blocks with an adaptor/probe arrangement:
% | ISI | ADAPTOR | GAP(?) | PROBE |
% User can specify adaptor/probe characteristics (spectrum, duration, etc.)
% As well as a range of probe levels and locations

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

F.Stot_ht = 0;
F.Sind_ = .025;
F.S_lt = F.X_lt + F.X_wd + 0.01*F.scrn_wd;
F.S_btm = F.X_btm;
F.S_wd=round(.45*F.scrn_wd); F.S_ht = round(.9*F.scrn_ht);

%Figure window
if(~exist1('H.searchfig') &...
      ~exist1('H.ablfig') &...
      ~exist1('H.itdfig') &...
      ~exist1('H.ildfig') &...
      ~exist1('H.ildfreqfig') &...
      ~exist1('H.AdaptedSpace_fig'))
   H.AdaptedSpace_fig = figure('Units','pixels',...
      'Position',[F.S_lt F.S_btm F.S_wd F.S_ht],...
      'Name','',...
      'NumberTitle','off');

% Title Bar
F.Stot_ht = F.Stot_ht + F.Sind_;
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 1-F.Stot_ht 1 F.Sind_],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontSize',12,...
   'FontWeight','bold',...
   'String', 'Spatial Adaptation Control');

F.Stot_ht = F.Stot_ht + 1.5*F.Sind_;

% increment test number
H.AdaptedSpace_inc_testnum = uicontrol('Style','pushbutton',...
   'Units','normalized',...
   'Position',[0.01 1-F.Stot_ht .48 F.Sind_],...
   'BackgroundColor','yellow',...
   'ForegroundColor','blue',...
   'FontWeight','bold',...
   'String', 'INCREMENT testnum',...
   'Callback','setInfo_AdaptedSpace');

% Engage TEST
H.AdaptedSpace_engage = uicontrol('Style','pushbutton',...
   'Units','normalized',...
   'Position',[.51 1-F.Stot_ht .48 F.Sind_],...
   'BackgroundColor','green',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'ENGAGE',...
   'Callback','Engage_AdaptedSpace');

F.Stot_ht = F.Stot_ht + F.Sind_;

%Record Data?
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0.01 1-F.Stot_ht .25 F.Sind_],...
   'String', 'Record Data?');
H.AdaptedSpace_recorddataflag = uicontrol('Style','checkbox',...
   'Units','normalized',...
   'Position',[.26 1-F.Stot_ht .05 F.Sind_],...
   'Callback','setInfo_AdaptedSpace');
% FN.data display
H.AdaptedSpace_data_FN = uicontrol('Style','text',...
   'Units','normalized',...
   'ForegroundColor','blue',...
   'Position',[.31 1-F.Stot_ht .17 F.Sind_],...
   'String', FN.data);

%Reset
H.AdaptedSpace_reset = uicontrol('Style','pushbutton',...
   'Units','normalized',...
   'Position',[.51 1-F.Stot_ht .48 F.Sind_],...
   'BackgroundColor','m',...
   'ForegroundColor','black',...
   'FontWeight','bold',...
   'String', 'Reset',...
   'Callback','Reset_AdaptedSpace');

F.Stot_ht = F.Stot_ht + F.Sind_;

%Pause TEST
H.AdaptedSpace_pause = uicontrol('Style','togglebutton',...
   'Units','normalized',...
   'Position',[.01 1-F.Stot_ht .48 F.Sind_],...
   'BackgroundColor',[.5 1 .5],...
   'ForegroundColor','black',...
   'FontWeight','bold',...
   'String', 'Pause');
%Exit
H.exitspace_FPR = uicontrol('Style','togglebutton',...
   'Units','normalized',...
   'Position',[.51 1-F.Stot_ht .48 F.Sind_],...
   'BackgroundColor',[1 1 1],...
   'ForegroundColor','yellow',...
   'FontWeight','bold',...
   'String', 'EXIT',...
   'Callback','Exit_space_FPR');
return


%Status bar
F.tot_ht = F.tot_ht + F.ind_ht *2;
H.space_FPR_status = uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Status: setting parameters');

%Remaining Reps and Trials
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 F.space_ht-F.tot_ht 100 F.ind_ht],...
   'String', 'Remaining Reps:');
H.space_FPR_remreps = uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[102 F.space_ht-F.tot_ht 40 F.ind_ht],...
   'String', get(H.space_FPR_numreps,'String'));
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[144 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', 'Trials:');
H.space_FPR_remtrials = uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[206 F.space_ht-F.tot_ht 40 F.ind_ht],...
   'String', '');

% ADAPTOR PARAMS
uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.01 1-FStot_ht 0.4 F.Sind_],...
    'BackgroundColor','green',...
    'ForegroundColor','black',...
    'FontSize',12,...
    'FontWeight','bold',...
    'HorizontalPosition','left',...
    'String','Adaptor Params:');

% TEST PROBE PARAMS
uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.51 1-FStot_ht 0.4 F.Sind_],...
    'BackgroundColor','green',...
    'ForegroundColor','black',...
    'FontSize',12,...
    'FontWeight','bold',...
    'String','Probe Params:');

% TRIAL SEQUENCE PARAMS
uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.51 1-FStot_ht 0.4 F.Sind_],...
    'BackgroundColor','green',...
    'ForegroundColor','black',...
    'FontSize',12,...
    'FontWeight','bold',...
    'String','Probe Params:');

return


%reminder
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'ForegroundColor',[0 1 1],...
   'FontSize',12,...
   'FontWeight','bold',...
   'String', 'HRTFs in DSPs: use *.eq HRIRs');

%Specify HRTF and Location files
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Specify Files');

%Fully-cued HRTF file
F.tot_ht = F.tot_ht + F.ind_ht;
H.spacefilepb = uicontrol('Parent',H.space_FPRfig,...
   'Style','pushbutton',...
   'Units','normalized',...
   'Position',[2 F.space_ht-F.tot_ht 195 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select Fully-cued Filter File',...
   'Callback','setInfo_space_FPR');
uicontrol('Parent',H.space_FPRfig,...
   'Style','text',...
   'Units','normalized',...
   'Position',[200 F.space_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'Fully-cued Test');
H.spaceflag = uicontrol('Parent',H.space_FPRfig,...
   'Style','checkbox',...
   'Units','normalized',...
   'Position',[315 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'Value',XStimParams.space_flag,...
   'Callback','setInfo_space_FPR');
F.tot_ht = F.tot_ht + F.ind_ht;
H.spacefilehdr = uicontrol('Parent',H.space_FPRfig,...
   'Style','text',...
   'Position',[0 F.space_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'Fully-cued HRTF File:');
H.spacefile = uicontrol('Parent',H.space_FPRfig,...
   'Style','text',...
   'Position',[120 F.space_ht-F.tot_ht F.space_wd-120 F.ind_ht],...
   'HorizontalAlignment','left',...
   'ForegroundColor','blue',...
   'String', FN.space_std);
uicontrol('Parent',H.space_FPRfig,...
   'Style','text',...
   'Units','normalized',...
   'Position',[200 F.space_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'ABL-alone Test');
H.space_ABLaloneflag = uicontrol('Parent',H.space_FPRfig,...
   'Style','checkbox',...
   'Units','normalized',...
   'Position',[315 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'Value',XStimParams.ABLalone_flag,...
   'Callback','setInfo_space_FPR');

%ILDAlone filter file
F.tot_ht = F.tot_ht + F.ind_ht;
H.space_ildalonepb = uicontrol('Parent',H.space_FPRfig,...
   'Style','pushbutton',...
   'Units','normalized',...
   'Position',[2 F.space_ht-F.tot_ht 195 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select ILD-alone Filter File',...
   'Callback','setinfo_space_FPR');
uicontrol('Parent',H.space_FPRfig,...
   'Style','text',...
   'Units','normalized',...
   'Position',[200 F.space_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'ILD-alone Test');
H.space_ildaloneflag = uicontrol('Parent',H.space_FPRfig,...
   'Style','checkbox',...
   'Units','normalized',...
   'Position',[315 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'Value',XStimParams.ildalone_flag,...
   'Callback','setInfo_space_FPR');
F.tot_ht = F.tot_ht + F.ind_ht;
H.space_ildalonefilehdr = uicontrol('Parent',H.space_FPRfig,...
   'Style','text',...
   'Position',[0 F.space_ht-F.tot_ht 100 F.ind_ht],...
   'String', 'ILD-alone File:');
H.space_ildalonefile = uicontrol('Parent',H.space_FPRfig,...
   'Style','text',...
   'Position',[100 F.space_ht-F.tot_ht F.space_wd-100 F.ind_ht],...
   'HorizontalAlignment','left',...
   'ForegroundColor','blue',...
   'String', FN.ildalone_std);

%ITDAlone filter file
F.tot_ht = F.tot_ht + F.ind_ht;
H.space_itdalonepb = uicontrol('Parent',H.space_FPRfig,...
   'Style','pushbutton',...
   'Units','normalized',...
   'Position',[2 F.space_ht-F.tot_ht 195 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select ITD-alone Filter File',...
   'Callback','setinfo_space_FPR');
uicontrol('Parent',H.space_FPRfig,...
   'Style','text',...
   'Units','normalized',...
   'Position',[200 F.space_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'ITD-alone Test');
H.space_itdaloneflag = uicontrol('Parent',H.space_FPRfig,...
   'Style','checkbox',...
   'Units','normalized',...
   'Position',[315 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'Value',XStimParams.itdalone_flag,...
   'Callback','setInfo_space_FPR');
F.tot_ht = F.tot_ht + F.ind_ht;
H.space_itdalonefilehdr = uicontrol('Parent',H.space_FPRfig,...
   'Style','text',...
   'Position',[0 F.space_ht-F.tot_ht 100 F.ind_ht],...
   'String', 'ITD-alone File:');
H.space_itdalonefile = uicontrol('Parent',H.space_FPRfig,...
   'Style','text',...
   'Position',[100 F.space_ht-F.tot_ht F.space_wd-100 F.ind_ht],...
   'HorizontalAlignment','left',...
   'ForegroundColor','blue',...
   'String', FN.itdalone_std);

F.tot_ht = F.tot_ht + F.ind_ht*2;
%Set Stimulus Parameters
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Set Stimulus Parameters');
F.tot_ht = F.tot_ht + F.ind_ht;    

% initialize stimFNs
H.space_FPR_initFN = uicontrol('Parent',H.space_FPRfig,...
   'Style','pushbutton',...
   'Units','normalized',...
   'Position',[2 F.space_ht-F.tot_ht 195 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Initialize stim FNs',...
   'Callback','setinfo_space_FPR');

H.space_FPR_DynamicRangepb = uicontrol('Style','togglebutton',...
   'Units','normalized',...
   'Position',[222 F.space_ht-F.tot_ht F.space_wd-222 F.ind_ht],...
   'BackgroundColor','yellow',...
   'ForegroundColor','blue',...
   'String', 'push for High Dynamic Range',...
   'Callback','setinfo_space_FPR');

F.tot_ht = F.tot_ht + F.ind_ht;
F.tot_ht = F.tot_ht + F.ind_ht;
F.tot_ht = F.tot_ht + F.ind_ht;
F.tot_ht = F.tot_ht + F.ind_ht;
% list stimFNs
for i = 1:4
H.space_FPR_stimFN(i) = uicontrol('Parent',H.space_FPRfig,...
   	'Style','text',...
   	'Units','normalized',...
   	'Position',[2+(i-1)*((F.space_wd-4)/4) F.space_ht - F.tot_ht (F.space_wd-4)/4 F.ind_ht*4],...
    'ForegroundColor','blue',...
   	'String', '');
end

%ABL (atten)
F.tot_ht = F.tot_ht + 2*F.ind_ht;
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'Atten (dB):');
H.space_FPR_ABL = uicontrol('Style','edit',...
   'Units','normalized',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_ABL),...
   'Callback','SetInfo_space_FPR');
% hold ITD?
H.space_FPR_ITDtxt = uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[165 F.space_ht-F.tot_ht 115 F.ind_ht],...
   'String', 'Hold ITD constant',...
   'BackgroundColor',[.8 .8 .8]);
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[300 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', 'ITD (us):');
H.space_FPR_ITD = uicontrol('Style','edit',...
   'Units','normalized',...
   'Position',[360 F.space_ht-F.tot_ht 40 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.curr_ITD),...
   'Callback','SetInfo_space_FPR');

%DUR
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'Stim Dur (ms):');
H.space_FPR_DUR = uicontrol('Style','edit',...
   'Units','normalized',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_stimdur),...
   'Callback','SetInfo_space_FPR');

%Hold ILD?
H.space_FPR_ILDtxt = uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[165 F.space_ht-F.tot_ht 115 F.ind_ht],...
   'String', 'Hold ILD constant',...
   'BackgroundColor',[.8 .8 .8]);
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[300 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'String', 'El');
H.space_FPR_el = uicontrol('Style','edit',...
   'Units','normalized',...
   'Position',[320 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.el),...
   'Callback','SetInfo_space_FPR');
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[350 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'String', 'Az');
H.space_FPR_az = uicontrol('Style','edit',...
   'Units','normalized',...
   'Position',[370 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.az),...
   'Callback','SetInfo_space_FPR');

%ISI
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'ISI (ms):');
H.space_FPR_ISI = uicontrol('Style','edit',...
   'Units','normalized',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.test_ISI));
%NumReps
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', '# Reps:');
H.space_FPR_numreps = uicontrol('Style','edit',...
   'Units','normalized',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.numreps),...
   'Callback','setinfo_space_FPR');

uicontrol('Style','text',...
   'Units','normalized',...
   'ForegroundColor','blue',...
   'BackgroundColor',[.8 .8 .8],...
   'Position',[150 F.space_ht-F.tot_ht 250 F.ind_ht],...
   'String', 'Plays each file, then repeats for nReps');

F.tot_ht = F.tot_ht + F.ind_ht*2;

%%%%%%%%%%%
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','normalized',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Run Test');

F.tot_ht = F.tot_ht + F.ind_ht;

end %end GUI specification

% Initialize Application and get AP2 and XBUS locks
if(S232('S2init', 0, 'INIT_PRIMARY', 20000) == 0)
   disp('Cannot initialize a primary process')
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

SetInfo_space_FPR;