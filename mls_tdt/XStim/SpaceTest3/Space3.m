%Space: the GUI for Space tests
% plays through DSPs that contain HRIRs for spatial locations


%*******************************************************************************
%	The Graphical User Interface
%*******************************************************************************

XStimParams.reset_flag = 0;
GUI.locations1 = [];
GUI.locations2 = [];
F.tot_ht = 0;
F.ind_ht = 20; 

pos_incr = 20;
%Figure window
if(~exist1('H.searchfig') &...
      ~exist1('H.ablfig') &...
      ~exist1('H.itdfig') &...
      ~exist1('H.ildfig') &...
      ~exist1('H.ildfreqfig') &...
      ~exist1('H.spacefig3'))
   Xstim_control_pos = get(H.fig,'Position');
   spacefig_l = Xstim_control_pos(3) - 300; F.space_wd = 400; F.space_ht = 950;
   H.spacefig3 = figure('Units','pixels',...
      'Position',[spacefig_l 0 F.space_wd F.space_ht],...
      'Name','Space3: DBL buffered virtual auditory space',...
      'NumberTitle','off');
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'ForegroundColor',[0 1 1],...
   'FontSize',12,...
   'FontWeight','bold',...
   'String', 'Stim written to disk before playout');

  %reminder
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'ForegroundColor',[0 1 1],...
   'FontSize',12,...
   'FontWeight','bold',...
   'String', 'HRIRs in DSPs: use *.eq HRIRs');

%Specify HRTF and Location files
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Specify Files');

%Space filter file (fully cued)
F.tot_ht = F.tot_ht + F.ind_ht;
H.spacefilepb = uicontrol('Parent',H.spacefig3,...
   'Style','pushbutton',...
   'Units','pixels',...
   'Position',[2 F.space_ht-F.tot_ht 195 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select Fully-Cued Filter File',...
   'Callback','setinfo_space3');
uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Units','pixels',...
   'Position',[200 F.space_ht-F.tot_ht 110 F.ind_ht],...
   'String', 'Fully-Cued Test');
H.spaceflag = uicontrol('Parent',H.spacefig3,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[315 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'Value',(XStimParams.space_flag),...
   'Callback','setinfo_space3');
F.tot_ht = F.tot_ht + F.ind_ht;
H.spacefilehdr = uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Position',[0 F.space_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'Fully-Cued File:');
H.spacefile = uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Position',[120 F.space_ht-F.tot_ht F.space_wd-120 F.ind_ht],...
   'HorizontalAlignment','left',...
   'ForegroundColor','blue',...
   'String', FN.space_eq);

%ILDAlone filter file
F.tot_ht = F.tot_ht + F.ind_ht;
H.ildalonepb = uicontrol('Parent',H.spacefig3,...
   'Style','pushbutton',...
   'Units','pixels',...
   'Position',[2 F.space_ht-F.tot_ht 195 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select ILDAlone Filter File',...
   'Callback','setinfo_space3');
uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Units','pixels',...
   'Position',[200 F.space_ht-F.tot_ht 110 F.ind_ht],...
   'String', 'ILDAlone Test');
H.ildaloneflag = uicontrol('Parent',H.spacefig3,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[315 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'Value',XStimParams.ildalone_flag,...
   'Callback','setinfo_space3');
F.tot_ht = F.tot_ht + F.ind_ht;
H.ildalonefilehdr = uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Position',[0 F.space_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'ILDAlone File:');
H.ildalonefile = uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Position',[120 F.space_ht-F.tot_ht F.space_wd-120 F.ind_ht],...
   'HorizontalAlignment','left',...
   'ForegroundColor','blue',...
   'String', FN.ildalone_eq);

%ITDAlone filter file
F.tot_ht = F.tot_ht + F.ind_ht;
H.itdalonepb = uicontrol('Parent',H.spacefig3,...
   'Style','pushbutton',...
   'Units','pixels',...
   'Position',[2 F.space_ht-F.tot_ht 195 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Select ITDAlone Filter File',...
   'Callback','setinfo_space3');
uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Units','pixels',...
   'Position',[200 F.space_ht-F.tot_ht 110 F.ind_ht],...
   'String', 'ITDAlone Test');
H.itdaloneflag = uicontrol('Parent',H.spacefig3,...
   'Style','checkbox',...
   'Units','pixels',...
   'Position',[315 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'Value',XStimParams.itdalone_flag,...
   'Callback','setInfo_space3');
F.tot_ht = F.tot_ht + F.ind_ht;
H.itdalonefilehdr = uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Position',[0 F.space_ht-F.tot_ht 100 F.ind_ht],...
   'String', 'ITDAlone File:');
H.itdalonefile = uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Position',[100 F.space_ht-F.tot_ht F.space_wd-100 F.ind_ht],...
   'HorizontalAlignment','left',...
   'ForegroundColor','blue',...
   'String', FN.itdalone_eq);

F.tot_ht = F.tot_ht + F.ind_ht *2;
% BLUE BAR: Set parameters
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Set Parameters');
F.tot_ht = F.tot_ht + F.ind_ht;

uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Units','pixels',...
   'Position',[2 F.space_ht-F.tot_ht 110 F.ind_ht],...
   	'BackgroundColor','red',...
   	'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Choose Stim Type');
H.space3_stim_type = uicontrol('Parent',H.spacefig3,...
   'Style','popup',...
   'Units','pixels',...
   'Position',[112 F.space_ht-F.tot_ht 110 F.ind_ht],...
   'BackgroundColor','white',...
   'String',...
   'Tone|Gammatone|Tone Stack|GammaTone Stack|BandPassed Noise|NarrowBand Noise|NarrowBand Noise Stack|BroadBand Noise|File|None',...
   'Value',get(H.stim_type,'Value'),...
   'Callback','setinfo_space3');

% allow re-choice of stimulus filename
H.space3_stimFNpb = uicontrol('Parent',H.spacefig3,...
   'Style','pushbutton',...
   'Units','pixels',...
   'Position',[224 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'FN',...
   'Callback','setinfo_space3');
H.space3_stimFN = uicontrol('Parent',H.spacefig3,...
   'Style','text',...
   'Units','pixels',...
   'ForegroundColor','blue',...
   'Position',[246 F.space_ht-F.tot_ht 154 F.ind_ht],...
   'String', FN.stim);


F.tot_ht = F.tot_ht + F.ind_ht*2;
% curr_freq
H.space3_freqtxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'center freq (Hz):');
H.space3_freq = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_freq),...
   'Callback','setInfo_space3');

F.tot_ht = F.tot_ht + F.ind_ht;
%ITD
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'ITD (us):');

H.space3_ITD = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_ITD),...
   'Callback','SetInfo_space3');

uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[165 F.space_ht-F.tot_ht 115 F.ind_ht],...
   'String', 'Hold ITD constant?');
H.ITDchk = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[280 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'Value', ~XStimParams.space_flag,...
   'Callback','SetInfo_space3');
%ABL
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'atten (dB):');

H.space3_ABL = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_ABL),...
   'Callback','SetInfo_space3');

%DUR
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'Stim Dur (ms):');
H.space3_DUR = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.curr_stimdur),...
   'Callback','SetInfo_space3');

%Hold ILD?
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[165 F.space_ht-F.tot_ht 115 F.ind_ht],...
   'String', 'Hold ILD constant?');
H.ILDchk = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[280 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'Value', XStimParams.itdalone_flag);
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[300 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'String', 'Az');
H.az = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[320 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.az),...
   'Callback','SetInfo_space3');
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[350 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'String', 'El');
H.el = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[370 F.space_ht-F.tot_ht 30 F.ind_ht],...
   'BackgroundColor','white',...
   'String', num2str(XStimParams.el),...
   'Callback','SetInfo_space3');

%ISI
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'ISI (ms):');

H.space3_ISI = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.test_ISI),...
   'Callback','SetInfo_space3');

%NumReps
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', '# Reps:');

H.space3_numreps = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[82 F.space_ht-F.tot_ht 60 F.ind_ht],...
   'String', num2str(XStimParams.numreps),...
   'Callback','SetInfo_space3');

   % modulate  sound
    F.tot_ht = F.tot_ht + F.ind_ht*2;
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
    uicontrol('Style','text',...
        'Units','pixels',...
        'FontWeight','bold',...
        'Position',[0 F.space_ht - F.tot_ht 80 F.ind_ht],...
        'String', '   Modulation:');
    H.Spacetest3mod_type = uicontrol('Style','popup',...
        'Units','pixels',...
        'Position',[80 F.space_ht-F.tot_ht 100 F.ind_ht],...
        'BackgroundColor','White',...
        'String',...
        'Tone|Sq wave|LP Noise|File|None',...
        'Value',mod_num,...
        'Callback','setInfo_Space3');
    H.Spacetest3mod_txt = uicontrol('Style','text',...
        'Units','pixels',...
        'HorizontalAlignment','right',...
        'Position',[180 F.space_ht-F.tot_ht 40 F.ind_ht],...
        'String', 'Freq:');
    H.Spacetest3mod_freq = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[220 F.space_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_freq(1)),...
        'Callback','setInfo_Space3');
    
    H.Spacetest3mod_txtA = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[260 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'depth');
    H.Spacetest3mod_depth = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[290 F.space_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_depth(1)),...
        'Callback','setInfo_Space3');
    F.tot_ht = F.tot_ht + F.ind_ht;

    H.Spacetest3mod_txtB = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[260 F.space_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'phase');
    H.Spacetest3mod_phase = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[290 F.space_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_phase(1)),...
        'Callback','setInfo_Space3');

    H.Spacetest3mod_pb = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[200 F.space_ht-F.tot_ht 25 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','blue',...
        'String', 'FN',...
        'Visible','off',...
        'Callback','setInfo_Space3');
 F.tot_ht = F.tot_ht + F.ind_ht;
   H.Spacetest3mod_txtC = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[10 F.space_ht-F.tot_ht 150 1.8*F.ind_ht],...
        'String', 'pi/2 phase begins at onset of squarewave max');

    % increment test number
F.tot_ht = F.tot_ht + F.ind_ht;
H.inc_testnum = uicontrol('Style','pushbutton',...
   'Units','pixels',...
   'Position',[250 F.space_ht-F.tot_ht 140 F.ind_ht],...
   'BackgroundColor','yellow',...
   'ForegroundColor','blue',...
   'FontWeight','bold',...
   'String', 'INCREMENT testnum',...
   'Callback','setInfo_space3');

% Only build files?? (do not play them out)
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[15 F.space_ht-F.tot_ht 120 F.ind_ht],...
   'String', 'ONLY Build Files?');
H.buildOnly = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[135 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'Callback','setInfo_space3');

% BLUE BAR: Run Test
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Run Test');

%Engage TEST
F.tot_ht = F.tot_ht + F.ind_ht;
H.engagespace = uicontrol('Style','pushbutton',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','green',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'ENGAGE',...
   'Callback','Engage_space3');


%Record Data?
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[85 F.space_ht-F.tot_ht 100 F.ind_ht],...
   'String', 'Record Data?');
H.recorddata = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[185 F.space_ht-F.tot_ht 20 F.ind_ht],...
   'Callback','setInfo_space3');

%Reset
H.resetspace = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[F.space_wd-100 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','cyan',...
   'ForegroundColor','black',...
   'FontWeight','bold',...
   'String', 'Reset',...
   'callback','Reset_space3');

% Data FN
F.tot_ht = F.tot_ht + F.ind_ht;
H.space3_recorddata_FN = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[85 F.space_ht-F.tot_ht 100 F.ind_ht],...
   'ForegroundColor','blue',...
   'String', FN.data);

%Pause TEST
F.tot_ht = F.tot_ht + F.ind_ht + pos_incr;
H.pausespace = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Pause');

%Exit
H.exitspace = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[F.space_wd-100 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'BackgroundColor','red',...
   'ForegroundColor','yellow',...
   'FontWeight','bold',...
   'String', 'EXIT',...
   'Callback','Exit_space3');

%Status bar
F.tot_ht = F.tot_ht + F.ind_ht *2;
H.status = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht F.space_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Status: setting parameters');

%Build/Play Trials
F.tot_ht = F.tot_ht + F.ind_ht + pos_incr;
H.buildplay = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 175 F.ind_ht],...
   'BackgroundColor','yellow',...
   'ForegroundColor','blue',...
   'String', 'Build/Play Status');

%Remaining Reps and Trials
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 100 F.ind_ht],...
   'String', 'Remaining reps:');
H.remreps = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[102 F.space_ht-F.tot_ht 40 F.ind_ht],...
   'String', get(H.space3_numreps,'String'));
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[144 F.space_ht-F.tot_ht 80 F.ind_ht],...
   'String', 'Trials:');
H.remtrials = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[226 F.space_ht-F.tot_ht 40 F.ind_ht],...
   'String', '');

%Plotting
F.tot_ht = F.tot_ht + F.ind_ht + pos_incr;
H.plotpsd = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[0 F.space_ht-F.tot_ht 175 F.ind_ht],...
   'BackgroundColor','yellow',...
   'ForegroundColor','black',...
   'FontWeight','bold',...
   'String', 'Plot PSD');

end %end GUI specification

stimuli_dir = FN.temp_stim_path;

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

SetInfo_space3;