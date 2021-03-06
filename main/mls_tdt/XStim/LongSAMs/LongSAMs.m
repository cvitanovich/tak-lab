%LongSAMs: the GUI for LongSAMs tests

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
XStimParams.LongSAMs_flag = 1;

% make tempStruct for play2_record2b_SAMB
tempStruct = struct('mdepth', XStimParams.mod_depth(1), ...
    'latten', abs(XStimParams.curr_ABL)-20, ...
    'ratten', abs(XStimParams.curr_ABL)-20, ...
    'SamplingPeriod', 1000000/TDT.Fs, ...
    'nptsTotalPlay', TDT.Fs, ...
    'buf_pts', TDT.Fs, ...
    'nCycles_perRep', 1, ...
    'ISI', XStimParams.test_ISI, ...
    'record', 0, ...
    'record_spikes', 1, ...
    'ramp_flag', 0, ...
    'noiScale', 10000);

F.tot_ht = 0;
F.ind_ht = 20; 
%Figure window
if(~exist1('H.searchfig') &...
        ~exist1('H.ablfig') &...
        ~exist1('H.itdfig') &...
        ~exist1('H.ildfig') &...
        ~exist1('H.ildfreqfig') &...
        ~exist1('H.LongSAMsfig'))
    temp = get(H.fig,'Position');
    F.l = temp(3) - 300; 
    F.wd = 400; 
    F.ht = 900;
    H.LongSAMsfig = figure('Units','pixels',...
        'Position',[F.l 100 F.wd F.ht],...
        'Name','Virtual Auditory  LongSAMs Measurement',...
        'NumberTitle','off');
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht F.wd F.ind_ht],...
        'ForegroundColor',[0 1 1],...
        'FontSize',12,...
        'FontWeight','bold',...
        'String', 'Stimuli played through dll.');
    
    %reminder
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht F.wd F.ind_ht],...
        'ForegroundColor',[0 1 1],...
        'FontSize',12,...
        'FontWeight','bold',...
        'String', 'HRTFs in DSPs: use *.eq HRIRs');
    
    %Specify HRTF and Location files
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht F.wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Specify Files');
    
    %Fully-cued HRTF file
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.LongSAMsfilepb = uicontrol('Parent',H.LongSAMsfig,...
        'Style','pushbutton',...
        'Units','pixels',...
        'Position',[2 F.ht-F.tot_ht 195 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Select Fully-cued Filter File',...
        'Callback','setInfo_LongSAMs');
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.LongSAMsfilehdr = uicontrol('Parent',H.LongSAMsfig,...
        'Style','text',...
        'Position',[0 F.ht-F.tot_ht 120 F.ind_ht],...
        'String', 'Fully-cued HRTF File:');
    H.LongSAMsfile = uicontrol('Parent',H.LongSAMsfig,...
        'Style','text',...
        'Position',[120 F.ht-F.tot_ht F.wd-120 F.ind_ht],...
        'HorizontalAlignment','left',...
        'ForegroundColor','blue',...
        'String', FN.space_eq);
    
    F.tot_ht = F.tot_ht + F.ind_ht*2;
    %Set Stimulus Parameters
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht F.wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Set Stimulus Parameters');
    F.tot_ht = F.tot_ht + F.ind_ht;    
    
    
    uicontrol('Parent',H.LongSAMsfig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[2 F.ht-F.tot_ht 210 F.ind_ht],...
        'BackgroundColor',[.8 .8 .8],...
        'ForegroundColor','black',...
        'FontWeight','bold',...
        'String', 'Stimulus is a BBnoise');
    
    F.tot_ht = F.tot_ht + F.ind_ht;
    if 0
    H.LongSAMs_DynamicRangepb = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[222 F.ht-F.tot_ht F.wd-222 F.ind_ht],...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'String', 'push for High Dynamic Range',...
        'Callback','setinfo_LongSAMs');
    
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.LongSAMs_DynamicRangeTxt = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[222 F.ht-F.tot_ht F.wd-222 F.ind_ht],...
        'BackgroundColor',[.8 .8 .8],...
        'String', 'SPL now 20dB below normal');
end

    %ABL (atten)
    F.tot_ht = F.tot_ht + 2*F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht 80 F.ind_ht],...
        'String', 'Atten (dB):');
    H.LongSAMs_ABL = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[82 F.ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.curr_ABL),...
        'Callback','SetInfo_LongSAMs');
    
    %DUR
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht 80 F.ind_ht],...
        'String', 'Stim Dur (ms):');
    H.LongSAMs_DUR = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[82 F.ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.curr_stimdur),...
        'Callback','SetInfo_LongSAMs');
    
    %ISI
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht 80 F.ind_ht],...
        'String', 'ISI (ms):');
    H.LongSAMs_ISI = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[82 F.ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.test_ISI));
    
    %NumReps
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht 80 F.ind_ht],...
        'String', '# Reps:');
    H.LongSAMs_numreps = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[82 F.ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.numreps),...
        'Callback','setinfo_LongSAMs');
    
    %%%%%%%%%%%%%%
        % modulate  sound
    F.tot_ht = F.tot_ht + F.ind_ht*2;
    
    if 0
        H.LongSAMstestmod_1or2 = uicontrol('Style','togglebutton',...
            'Units','pixels',...
            'Position',[0 F.ht-F.tot_ht 160 F.ind_ht],...
            'BackgroundColor','yellow',...
            'ForegroundColor','blue',...
            'FontWeight','bold',...
            'String', 'push for single mod',...
            'value',0,...
            'Callback','setinfo_LongSAMs');
    end
    
    H.LongSAMstestmod_1or2 = uicontrol('Style','popup',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'value', 1,...
        'String', ['Double|Single|Preset 4 freqs|Preset 7 freqs'],...
        'Callback','setinfo_LongSAMs');

        
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.LongSAMstestmod_txt(1) = uicontrol('Style','text',...
        'Units','pixels',...
        'FontWeight','bold',...
        'Position',[0 F.ht - F.tot_ht 150 F.ind_ht],...
        'String', '  Dbl Modulated- fast:');
    H.LongSAMstestmod_txt(4) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'HorizontalAlignment','right',...
        'Position',[150 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', 'Freq:');
    H.LongSAMstestmod_freq(1) = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[190 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_freq(1)),...
        'Callback','setInfo_LongSAMs');
        
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'FontWeight','bold',...
        'Position',[0 F.ht - F.tot_ht 280 F.ind_ht],...
        'String', 'Set period radix 2, <= 100 & > slow freq');
    
    uicontrol('Style','text',...
        'Units','pixels',...
        'BackgroundColor',[.8 .8 .8],...
        'ForegroundColor',[.6 0 .6],...
         'FontSize',12,...
      'Position',[300 F.ht - F.tot_ht 100 F.ind_ht*3.5],...
        'String', 'for 4x,7x mods set lowest freq (1/16Hz?) from single freq');
    
%%%%%%%%%%%%%%%%%%%
    % modulate  sound
    F.tot_ht = F.tot_ht + F.ind_ht;
    
    uicontrol('Style','text',...
        'Units','pixels',...
        'FontWeight','bold',...
        'Position',[0 F.ht - F.tot_ht 150 F.ind_ht],...
        'String', '  Modulated slowly');
    uicontrol('Style','text',...
        'Units','pixels',...
        'HorizontalAlignment','right',...
        'Position',[150 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', 'Freq:');
    % linlog indicated by mod_type
    H.LongSAMstestmod_linlog(1) = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[150 F.ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'FontWeight','bold',...
        'String', 'linear',...
        'value',0,...
        'Callback','setinfo_LongSAMs');
    H.LongSAMstestmod_txt(6) = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[230 F.ht-F.tot_ht 30 F.ind_ht],...
        'String', 'dB');
    H.LongSAMstestmod_dB(1) = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[260 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_dB(1)),...
        'Callback','setInfo_LongSAMs');

    F.tot_ht = F.tot_ht + F.ind_ht;    
    H.LongSAMstestmod_txt(7) = uicontrol(        'Style','text',...
        'Units','pixels',...
        'HorizontalAlignment','right',...
        'Position',[150 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', 'Freq:');
    H.LongSAMstestmod_freq(2) = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[190 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_freq(2)),...
        'Callback','setInfo_LongSAMs');
    
    H.LongSAMstestmod_txt(8) = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[230 F.ht-F.tot_ht 30 F.ind_ht],...
        'String', 'depth');
    H.LongSAMstestmod_depth(1) = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[260 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_depth(1)),...
        'Callback','setInfo_LongSAMs');
    
    H.LongSAMstestmod_txt(9) = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[300 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', 'nCycles');
    
    H.LongSAMstestmod_nCycles_perRep = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[340 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(tempStruct.nCycles_perRep),...
        'Callback','setInfo_LongSAMs');
    
    
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'FontWeight','bold',...
        'Position',[0 F.ht - F.tot_ht 280 F.ind_ht],...
        'String', 'Set period between .01 and 16 seconds');
    
    %%%%%%%%%%%
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht F.wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Run Test');
    
    F.tot_ht = F.tot_ht + F.ind_ht;
    
    % increment test number
    H.LongSAMs_inc_testnum = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[250 F.ht-F.tot_ht 140 F.ind_ht],...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'FontWeight','bold',...
        'String', 'INCREMENT testnum',...
        'Callback','setInfo_LongSAMs');
    
    %Engage TEST
    H.engageLongSAMs = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','green',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'ENGAGE',...
        'Callback','Engage_LongSAMs(tempStruct)');
    
    %Record Data?
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[85 F.ht-F.tot_ht 100 F.ind_ht],...
        'String', 'Record Data?');
    H.LongSAMs_recorddata = uicontrol('Style','checkbox',...
        'Units','pixels',...
        'Value', 1,...
        'Position',[185 F.ht-F.tot_ht 20 F.ind_ht],...
        'Callback','setInfo_LongSAMs');
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.LongSAMsdata_FN = uicontrol('Style','text',...
        'Units','pixels',...
        'ForegroundColor','blue',...
        'Position',[85 F.ht-F.tot_ht 100 F.ind_ht],...
        'String', FN.data);
    
    %Reset
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.resetLongSAMs = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[F.wd-100 F.ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','cyan',...
        'ForegroundColor','black',...
        'FontWeight','bold',...
        'String', 'Reset',...
        'Callback','Reset_LongSAMs');
    
    %Pause TEST
    F.ind_ht = 20; 
    F.tot_ht = F.tot_ht + F.ind_ht *2;
    H.pauseLongSAMs = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht 80 F.ind_ht*2],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Pause');
    
    %Exit
    H.exitLongSAMs = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[F.wd-100 F.ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','yellow',...
        'FontWeight','bold',...
        'String', 'EXIT',...
        'Callback','Exit_LongSAMs');
    
    %Status bar
    F.tot_ht = F.tot_ht + F.ind_ht *2;
    H.LongSAMs_status = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht F.wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Status: setting parameters');
    
    %Remaining Reps and Trials
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.ht-F.tot_ht 100 F.ind_ht],...
        'String', 'Remaining Reps:');
    H.LongSAMs_remreps = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[102 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', get(H.LongSAMs_numreps,'String'));
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[144 F.ht-F.tot_ht 60 F.ind_ht],...
        'String', 'Trials:');
    H.LongSAMs_remtrials = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[206 F.ht-F.tot_ht 40 F.ind_ht],...
        'String', '');
    
end %end GUI specification

SetInfo_LongSAMs;