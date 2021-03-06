%PE_env: the GUI for Precedence Effect with gammatoned envelopes tests
    % presents two stimuli, sets delay, sets gammatone from which to take
    % envelope, etc.

%*******************************************************************************
%	The Graphical User Interface
%*******************************************************************************

GUI.locations1 = [];
GUI.locations2 = [];
XStimParams.reset_flag = 0;
XStimParams.test_ISI = 1000;

XStimParams.space_flag == 0;          % for fully cued
XStimParams.delay = [1 2 3 5 10];
XStimParams.cF = [500 1000 2000];


global FigParams
FigParams.handle = '';

F.pos_incr = 20;
F.tot_ht = 0;
F.ind_ht = 20; 
%Figure window
if(~exist1('H.searchfig') &...
        ~exist1('H.ablfig') &...
        ~exist1('H.itdfig') &...
        ~exist1('H.ildfig') &...
        ~exist1('H.ildfreqfig') &...
        ~exist1('H.PE_envfig'))
    temp = get(H.fig,'Position');
    F.PE_env_l = temp(3) - 300; 
    F.PE_env_wd = 400; 
    F.PE_env_ht = 900;
    H.PE_envfig = figure('Units','pixels',...
        'Position',[F.PE_env_l 100 F.PE_env_wd F.PE_env_ht],...
        'Name','Precedence Effect with gammatoned envelopes',...
        'NumberTitle','off');
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht F.PE_env_wd F.ind_ht],...
        'ForegroundColor','blue',...
        'FontWeight','bold',...
        'String', 'PE_env');
    %reminder
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht F.PE_env_wd F.ind_ht],...
        'ForegroundColor','blue',...
        'String', 'spatial HRTFs loaded to DSPs with *.eq');
    
    %Specify HRTF and Location files
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht F.PE_env_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Specify Files');
    
    %Fully-cued HRTF file
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.spacefilepb = uicontrol('Parent',H.PE_envfig,...
        'Style','pushbutton',...
        'Position',[2 F.PE_env_ht-F.tot_ht 195 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Select Fully-cued Filter File',...
        'Callback','setInfo_PE_env');
    
    H.spaceflag = uicontrol('Parent',H.PE_envfig,...
        'Style','togglebutton',...
        'Position',[240 F.PE_env_ht-F.tot_ht 120 F.ind_ht],...
        'Value',(XStimParams.space_flag),...
        'String','Push for ABLequal',...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'Callback','setInfo_PE_env');

   F.tot_ht = F.tot_ht + F.ind_ht;

   H.spacefilehdr = uicontrol('Parent',H.PE_envfig,...
        'Style','text',...
        'Position',[0 F.PE_env_ht-F.tot_ht 120 F.ind_ht],...
        'String', 'Fully-cued HRTF File:');
    H.spacefile = uicontrol('Parent',H.PE_envfig,...
        'Style','text',...
        'Position',[120 F.PE_env_ht-F.tot_ht F.PE_env_wd-120 F.ind_ht],...
        'HorizontalAlignment','left',...
        'ForegroundColor','blue',...
        'String', FN.space_std);
            
    F.tot_ht = F.tot_ht + F.ind_ht *1.5;
    %Other Stimulus Parameters
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht F.PE_env_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Set Stimulus Parameters');
    F.tot_ht = F.tot_ht + F.ind_ht;
    
     
    % offset source2 from source1
    F.tot_ht = F.tot_ht + F.ind_ht*1.5;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht 155 F.ind_ht],...
        'FontWeight','bold',...
        'FontSize',10,...
        'Backgroundcolor',[.8 .8 .8], ...
        'String', 'loc2 = loc1 + OFFSET:');
  
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[155 F.PE_env_ht-F.tot_ht 25 F.ind_ht],...
        'String', 'elev');
    H.offset_el = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[180 F.PE_env_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.offset_el(2)),...
        'Callback','setinfo_PE_env');
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[220 F.PE_env_ht-F.tot_ht 25 F.ind_ht],...
        'String', 'azim');
    H.offset_az = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[250 F.PE_env_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.offset_az(2)),...
        'Callback','setinfo_PE_env');
    
 
    

    F.tot_ht = F.tot_ht + F.ind_ht;
        uicontrol('Style','text',...
        'Units','pixels',...
        'Backgroundcolor',[.8 .8 .8], ...
        'ForeGroundColor','blue', ...
        'Position',[270 F.PE_env_ht-F.tot_ht 130 F.ind_ht],...
        'String', 'silence periods (ms):');

    %ABL
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht 100 F.ind_ht],...
        'String', 'atten (dB)');
    H.ABL = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[100 F.PE_env_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.curr_ABL),...
        'Callback', 'setinfo_PE_env');
    
    if 0
    %ISI
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[180 F.PE_env_ht-F.tot_ht 60 F.ind_ht],...
        'String', 'ISI (ms)');
    H.ISI = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[240 F.PE_env_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.test_ISI),...
        'Callback', 'setinfo_PE_env');
end

    % silence_lead
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[290 F.PE_env_ht-F.tot_ht 70 F.ind_ht],...
        'String', 'leading');
    H.PE_env_silence_lead = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[360 F.PE_env_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.silence_lead),...
        'Callback', 'setinfo_PE_env');

    %DUR
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht 100 F.ind_ht],...
        'String', 'total stim dur (ms)');
    H.DUR = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[100 F.PE_env_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.curr_stimdur),...
        'Callback','setinfo_PE_env');
    
    %NumReps
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[180 F.PE_env_ht-F.tot_ht 60 F.ind_ht],...
        'String', '# Reps');
    H.numreps = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[240 F.PE_env_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.numreps),...
        'Callback','setinfo_PE_env');
    
        % silence_lead
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[290 F.PE_env_ht-F.tot_ht 70 F.ind_ht],...
        'String', 'trailing');
    H.PE_env_silence_trail = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[360 F.PE_env_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.silence_trail),...
        'Callback', 'setinfo_PE_env');


    % delays
    F.tot_ht = F.tot_ht + F.ind_ht*2;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht 100 F.ind_ht],...
        'String', 'delays (ms)');
    H.PE_env_delay = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[100 F.PE_env_ht-F.tot_ht 160 F.ind_ht],...
        'String', num2str(XStimParams.delay),...
        'Callback', 'setinfo_PE_env');
    
    % gammatone cFs
    F.tot_ht = F.tot_ht + F.ind_ht*2;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht 100 F.ind_ht],...
        'String', 'gammatone cF (Hz)');
    H.PE_env_cF = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[100 F.PE_env_ht-F.tot_ht 160 F.ind_ht],...
        'String', num2str(XStimParams.cF),...
        'Callback', 'setinfo_PE_env');
    
    F.tot_ht = F.tot_ht + F.ind_ht*6;

%%%%%%%%%%%%%% run control
    F.tot_ht = F.tot_ht + F.ind_ht*1.5;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht F.PE_env_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Run Control');
    
    %Engage TEST
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.engagePE_env = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','green',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'ENGAGE',...
        'Callback','Engage_PE_env');
    
    %Record Data?
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[85 F.PE_env_ht-F.tot_ht 100 F.ind_ht],...
        'String', 'Record Data?');
    H.recorddata = uicontrol('Style','checkbox',...
        'Units','pixels',...
        'Position',[185 F.PE_env_ht-F.tot_ht 20 F.ind_ht],...
        'Callback', 'setinfo_PE_env');
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.recorddata_FN = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[85 F.PE_env_ht-F.tot_ht 100 F.ind_ht],...
        'String', FN.data);
    
    % increment test number
    H.inc_testnum = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[250 F.PE_env_ht-F.tot_ht 140 F.ind_ht],...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'FontWeight','bold',...
        'String', 'INCREMENT testnum',...
        'Callback','setInfo_PE_env');
    
    %Pause TEST
    F.ind_ht = 20; 
    F.tot_ht = F.tot_ht + F.ind_ht + F.pos_incr;
    H.pausePE_env = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Pause');
   
    %Reset
    H.resetPE_env = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[200 F.PE_env_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','cyan',...
        'ForegroundColor','black',...
        'FontWeight','bold',...
        'String', 'Reset',...
        'Callback','Reset_PE_env');
    
    %Exit
    H.exitPE_env = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[F.PE_env_wd-100 F.PE_env_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','yellow',...
        'FontWeight','bold',...
        'String', 'EXIT',...
        'Callback','Exit_PE_env');
    
    %Build/Play Trials
    F.ind_ht = 20; 
    F.tot_ht = F.tot_ht + F.ind_ht + F.pos_incr;
    H.status = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht F.PE_env_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Status: setting parameters');
    
    %Remaining Reps and Trials
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.PE_env_ht-F.tot_ht 120 F.ind_ht],...
        'String', 'Remaining reps:');
    H.remreps = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[122 F.PE_env_ht-F.tot_ht 40 F.ind_ht],...
        'String', get(H.numreps,'String'));
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[164 F.PE_env_ht-F.tot_ht 80 F.ind_ht],...
        'String', 'Trials:');
    H.remtrials = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[246 F.PE_env_ht-F.tot_ht 40 F.ind_ht],...
        'String', '');
    
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

SetInfo_PE_env;