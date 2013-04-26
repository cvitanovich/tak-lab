%Mc2Source: the GUI for Mc2Source tests
    % presents two Dean & McAlpine-like stimuli

%*******************************************************************************
%	The Graphical User Interface
%*******************************************************************************

GUI.locations1 = [];
GUI.locations2 = [];
XStimParams.reset_flag = 0;
XStimParams.bandwidth = [0.3333 0.3333 0.3333];
XStimParams.bandlimits = [2000 11000 2000 11000 2000 11000];
XStimParams.test_ISI = 1000;

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
        ~exist1('H.Mc2Sourcefig'))
    temp = get(H.fig,'Position');
    F.Mc2Source_l = temp(3) - 300; 
    F.Mc2Source_wd = 400; 
    F.Mc2Source_ht = 900;
    H.Mc2Sourcefig = figure('Units','pixels',...
        'Position',[F.Mc2Source_l 100 F.Mc2Source_wd F.Mc2Source_ht],...
        'Name','TWO-source Space Measurement with Dean&McAlpine like stimuli',...
        'NumberTitle','off');
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd F.ind_ht],...
        'ForegroundColor','blue',...
        'FontWeight','bold',...
        'String', 'Mc TWO-sources');
    %reminder
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd F.ind_ht],...
        'ForegroundColor','blue',...
        'String', 'spatial HRTFs loaded to DSPs with *.eq');
    
    %Specify HRTF and Location files
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Specify Files');
    
    %Fully-cued or ABL-equal FC HRTF file
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.spacefilepb = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','pushbutton',...
        'Position',[2 F.Mc2Source_ht-F.tot_ht 195 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Select Fully-cued Filter File',...
        'Callback','setInfo_Mc2Source');
    
    H.spaceflag = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','togglebutton',...
        'Position',[240 F.Mc2Source_ht-F.tot_ht 120 F.ind_ht],...
        'Value',(XStimParams.space_flag),...
        'String','Push for ABLequal',...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'Callback','setInfo_Mc2Source');
    F.tot_ht = F.tot_ht + F.ind_ht;

   H.spacefilehdr = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht 120 F.ind_ht],...
        'String', 'Fully-cued HRTF File:');
    H.spacefile = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Position',[120 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd-120 F.ind_ht],...
        'HorizontalAlignment','left',...
        'ForegroundColor','blue',...
        'String', FN.space_std);
            
    F.tot_ht = F.tot_ht + F.ind_ht *1.5;
    %Other Stimulus Parameters
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Set Stimulus Parameters');
    F.tot_ht = F.tot_ht + F.ind_ht;
    
    %First Stimulus types: popupmenu
    uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht - F.tot_ht 60 F.ind_ht],...
        'ForegroundColor','blue',...
        'FontWeight','bold',...
        'String', '1st Stim:');
    H.Mc2Src_stim_type = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','popup',...
        'Units','pixels',...
        'Position',[60 F.Mc2Source_ht-F.tot_ht 120 F.ind_ht],...
        'BackgroundColor','white',...
        'String',...
        'Tone|Gammatone|Tone Stack|GammaTone Stack|BandPassed Noise|NarrowBand Noise|NarrowBand Noise Stack|BroadBand Noise|File|BBN-LNN|None',...
        'Value',get(H.stim_type,'Value'),...
        'Callback','setinfo_Mc2Source');
    
    % only one will be displayed, either
    H.Mc2Src_stim_filenamepb = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','pushbutton',...
        'Units','pixels',...
        'Position',[180 F.Mc2Source_ht-F.tot_ht 25 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','blue',...
        'String', 'FN',...
        'Callback','setinfo_Mc2Source');
    H.Mc2Src_stim_filename = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[205 F.Mc2Source_ht-F.tot_ht 105 F.ind_ht],...
        'ForegroundColor','blue',...
        'String', FN.stim);
    % or
    H.Mc2Src_freq_txt = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[180 F.Mc2Source_ht-F.tot_ht 90 F.ind_ht],...
        'Visible','off',...
        'String', 'Frequency');
    H.Mc2Src_freq = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','edit',...
        'Units','pixels',...
        'Visible','off',...
        'Position',[270 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.freq(1)),...
        'Callback','setinfo_Mc2Source');
    
    uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[340 F.Mc2Source_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'factor');
    H.Mc2Src_factor = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','edit',...
        'Units','pixels',...
        'Position',[370 F.Mc2Source_ht-F.tot_ht 30 F.ind_ht],...
        'String', num2str(XStimParams.factor(1)),...
        'Callback','setinfo_Mc2Source');
    
    % modulate 1st sound
    F.tot_ht = F.tot_ht + F.ind_ht;
    switch XStimParams.mod_type
        case 'Tone'
            mod_num = 1;
        case 'File'
            mod_num = 2;
        case 'LP Noise'
            mod_num = 3;
        case 'None'
            mod_num = 4;
        otherwise
            mod_num = 1;
    end
    uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'FontWeight','bold',...
        'Position',[0 F.Mc2Source_ht - F.tot_ht 80 F.ind_ht],...
        'String', '   Modulation:');
    H.Mc2Src_mod_type = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','popup',...
        'Units','pixels',...
        'Position',[80 F.Mc2Source_ht-F.tot_ht 100 F.ind_ht],...
        'BackgroundColor','White',...
        'String',...
        'Tone|LP Noise|File|None',...
        'Value',mod_num,...
        'Callback','setinfo_Mc2Source');
    H.Mc2Src_mod_txt = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'HorizontalAlignment','right',...
        'Position',[180 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
        'String', 'Mod Freq: ');
    H.Mc2Src_mod_freq = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','edit',...
        'Units','pixels',...
        'Position',[230 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_freq(1)),...
        'Callback','setinfo_Mc2Source');
    
    H.Mc2Src_mod_txtA = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[270 F.Mc2Source_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'depth');
    H.Mc2Src_mod_depth = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','edit',...
        'Units','pixels',...
        'Position',[300 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_depth(1)),...
        'Callback','setinfo_Mc2Source');
    H.Mc2Src_mod_txtB = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[340 F.Mc2Source_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'phase');
    H.Mc2Src_mod_phase = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','edit',...
        'Units','pixels',...
        'Position',[370 F.Mc2Source_ht-F.tot_ht 30 F.ind_ht],...
        'String', num2str(XStimParams.mod_phase(1)),...
        'Callback','setinfo_Mc2Source');
    
    H.Mc2Src_mod_pb = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','pushbutton',...
        'Units','pixels',...
        'Position',[200 F.Mc2Source_ht-F.tot_ht 25 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','blue',...
        'String', 'FN',...
        'Visible','off',...
        'Callback','setinfo_Mc2Source');
    
    %Second Stimulus types: popupmenu
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht - F.tot_ht 60 F.ind_ht],...
        'ForegroundColor','blue',...
        'FontWeight','bold',...
        'String', '2nd Stim:');
    H.Mc2Src_stim_type2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','popup',...
        'Units','pixels',...
        'Position',[60 F.Mc2Source_ht-F.tot_ht 120 F.ind_ht],...
        'BackgroundColor','white',...
        'String',...
        'Tone|Gammatone|Tone Stack|GammaTone Stack|BandPassed Noise|NarrowBand Noise|NarrowBand Noise Stack|BroadBand Noise|File|BBN-LNN|None',...
        'Value',get(H.stim_type2,'Value'),...
        'Callback','setinfo_Mc2Source');
    
    % only one will be visible, either:
    H.Mc2Src_stim_filename2pb = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','pushbutton',...
        'Units','pixels',...
        'Position',[180 F.Mc2Source_ht-F.tot_ht 25 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','blue',...
        'String', 'FN',...
        'Callback','setinfo_Mc2Source');
    H.Mc2Src_stim_filename2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[205 F.Mc2Source_ht-F.tot_ht 105 F.ind_ht],...
        'ForegroundColor','blue',...
        'String', FN.stim2);
    % or
    H.Mc2Src_freq2_txt = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[180 F.Mc2Source_ht-F.tot_ht 90 F.ind_ht],...
        'Visible','off',...
        'String', 'Frequency');
    H.Mc2Src_freq2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','edit',...
        'Units','pixels',...
        'Visible','off',...
        'Position',[270 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.freq(2)),...
        'Callback','setinfo_Mc2Source');
    
    uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[340 F.Mc2Source_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'factor');
    H.Mc2Src_factor2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','edit',...
        'Units','pixels',...
        'Position',[370 F.Mc2Source_ht-F.tot_ht 30 F.ind_ht],...
        'String', num2str(XStimParams.factor(2)),...
        'Callback','setinfo_Mc2Source');
    
    % modulate 2nd sound
    F.tot_ht = F.tot_ht + F.ind_ht;
    switch XStimParams.mod_type2
        case 'Tone'
            mod_num = 1;
        case 'File'
            mod_num = 2;
        case 'LP Noise'
            mod_num = 3;
        case 'None'
            mod_num = 4;
        otherwise
            mod_num = 1;
    end
    uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'FontWeight','bold',...
        'Position',[0 F.Mc2Source_ht - F.tot_ht 80 F.ind_ht],...
        'String', '  Modulation:');
    H.Mc2Src_mod_type2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','popup',...
        'Units','pixels',...
        'Position',[80 F.Mc2Source_ht-F.tot_ht 100 F.ind_ht],...
        'BackgroundColor','White',...
        'String',...
        'Tone|LP Noise|File|None',...
        'Value',mod_num,...
        'Callback','setinfo_Mc2Source');
    H.Mc2Src_mod_txt2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'HorizontalAlignment','right',...
        'Position',[180 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
        'String', 'Mod Freq:');
    H.Mc2Src_mod_freq2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','edit',...
        'Units','pixels',...
        'Position',[230 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_freq(2)),...
        'Callback','setinfo_Mc2Source');
    H.Mc2Src_mod_txtA2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[270 F.Mc2Source_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'depth');
    H.Mc2Src_mod_depth2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','edit',...
        'Units','pixels',...
        'Position',[300 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.mod_depth(2)),...
        'Callback','setinfo_Mc2Source');
    H.Mc2Src_mod_txtB2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','text',...
        'Units','pixels',...
        'Position',[340 F.Mc2Source_ht-F.tot_ht 30 F.ind_ht],...
        'String', 'phase');
    H.Mc2Src_mod_phase2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','edit',...
        'Units','pixels',...
        'Position',[370 F.Mc2Source_ht-F.tot_ht 30 F.ind_ht],...
        'String', num2str(XStimParams.mod_phase(2)),...
        'Callback','setinfo_Mc2Source');
    
    H.Mc2Src_mod_pb2 = uicontrol('Parent',H.Mc2Sourcefig,...
        'Style','pushbutton',...
        'Units','pixels',...
        'Position',[200 F.Mc2Source_ht-F.tot_ht 25 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','blue',...
        'String', 'FN',...
        'Visible','off',...
        'Callback','setinfo_Mc2Source');
    

    F.tot_ht = F.tot_ht + F.ind_ht;
    H.Mc2Src_DynamicRangepb = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[180 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd-180 F.ind_ht],...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'String', 'push for High Dynamic Range',...
        'Callback','setinfo_Mc2Source');
    
    H.Mc2Src_DynamicRangeTxt = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht 180 F.ind_ht],...
        'BackgroundColor',[.8 .8 .8],...
        'String', 'SPL now 20dB below normal');
   
    % offset source2 from source1
    F.tot_ht = F.tot_ht + F.ind_ht*1.5;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht 155 F.ind_ht],...
        'FontWeight','bold',...
        'FontSize',10,...
        'Backgroundcolor',[.8 .8 .8], ...
        'String', 'loc2 = loc1 + OFFSET:');
    
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[155 F.Mc2Source_ht-F.tot_ht 25 F.ind_ht],...
        'String', 'elev');
    H.offset_el = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[180 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.offset_el(2)),...
        'Callback','setinfo_Mc2Source');
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[220 F.Mc2Source_ht-F.tot_ht 25 F.ind_ht],...
        'String', 'azim');
    H.offset_az = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[250 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.offset_az(2)),...
        'Callback','setinfo_Mc2Source');
    
    % bandlimits
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.bandlimits_txt1 = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht 80 F.ind_ht],...
        'String', 'bandlimits  s1:');
    H.bandlimits1 = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[82 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
        'String', num2str(XStimParams.bandlimits(1)),...
        'Callback','setinfo_Mc2Source');
    H.bandlimits2 = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[132 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.bandlimits(2)),...
        'Callback','setinfo_Mc2Source');
    
    H.bandlimits_txt2 = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[174 F.Mc2Source_ht-F.tot_ht 24 F.ind_ht],...
        'String', 's2:');
    H.bandlimits3 = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[198 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
        'String', num2str(XStimParams.bandlimits(3)),...
        'Callback','setinfo_Mc2Source');
    H.bandlimits4 = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[248 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
        'String', num2str(XStimParams.bandlimits(4)),...
        'Callback','setinfo_Mc2Source');
    

    F.tot_ht = F.tot_ht + F.ind_ht;
        uicontrol('Style','text',...
        'Units','pixels',...
        'Backgroundcolor',[.8 .8 .8], ...
        'ForeGroundColor','blue', ...
        'Position',[270 F.Mc2Source_ht-F.tot_ht 130 F.ind_ht],...
        'String', 'silence periods (ms):');

    %ABL
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht 100 F.ind_ht],...
        'String', 'atten (dB)');
    H.ABL = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[100 F.Mc2Source_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.curr_ABL),...
        'Callback', 'setinfo_Mc2Source');
    
    %ISI
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[180 F.Mc2Source_ht-F.tot_ht 60 F.ind_ht],...
        'String', 'ISI (ms)');
    H.ISI = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[240 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.test_ISI),...
        'Callback', 'setinfo_Mc2Source');
    
    % silence_lead
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[290 F.Mc2Source_ht-F.tot_ht 70 F.ind_ht],...
        'String', 'leading');
    H.Mc2Source_silence_lead = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[360 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.silence_lead),...
        'Callback', 'setinfo_Mc2Source');

    %DUR
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht 100 F.ind_ht],...
        'String', 'total stim dur (ms)');
    H.DUR = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[100 F.Mc2Source_ht-F.tot_ht 60 F.ind_ht],...
        'String', num2str(XStimParams.curr_stimdur),...
        'Callback','setinfo_Mc2Source');
    
    %NumReps
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[180 F.Mc2Source_ht-F.tot_ht 60 F.ind_ht],...
        'String', '# Reps');
    H.numreps = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[240 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.numreps),...
        'Callback','setinfo_Mc2Source');
    
        % silence_lead
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[290 F.Mc2Source_ht-F.tot_ht 70 F.ind_ht],...
        'String', 'trailing');
    H.Mc2Source_silence_trail = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[360 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', num2str(XStimParams.silence_trail),...
        'Callback', 'setinfo_Mc2Source');

%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% SETup SPLs  
F.tot_ht = F.tot_ht + F.ind_ht*1.5;
H.Mc2Src_txt(20) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd F.ind_ht],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'SetUp SPLs');
  
    
    
    %%%%%%%%%% SRC2 FOCAL SPLs
F.tot_ht = F.tot_ht + F.ind_ht;
H.Mc2Src_txt(18) = uicontrol('Style','text',...
   'Units','pixels',...
   'BackgroundColor',[.8 .8 .8],...
   'ForegroundColor','blue',...
   'FontSize',10,...
   'Position',[0 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd-70 F.ind_ht],...
   'String', 'S1: constant probe. Set s2 (adapter) focal-range SPL:');
H.Mc2Src_txt(19) = uicontrol('Style','togglebutton',...
    'Units','pixels',...
    'Position',[F.Mc2Source_wd-70 F.Mc2Source_ht-F.tot_ht 70 F.ind_ht],...
    'BackgroundColor','yellow',...
    'ForegroundColor','blue',...
    'String', 'swap S1/S2',...
    'Callback','setinfo_Mc2Source');

%SPL FOCAL-range lower-SPL slider   (slider#1)
F.tot_ht = F.tot_ht + F.ind_ht;
H.Mc2Source_SPLslider(1) = uicontrol('Style','slider',...
   'Units','pixels',...
   'Position',[0 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd-60 F.ind_ht *.75],...
   'BackgroundColor','blue',...
   'Min',0,'Max',30,...
   'SliderStep',[1/31 5/31],...
   'Value',XStimParams.focalSPL(1),...
   'Callback','setinfo_Mc2Source');
H.Mc2Source_SPLslider_txt(1) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[F.Mc2Source_wd-60 F.Mc2Source_ht-F.tot_ht 60 F.ind_ht*.75],...
   'FontSize',8,...
   'String', [num2str(XStimParams.focalSPL(1)) ' dB']);
%SPL FOCAL-range upper-SPL slider  (slider#2)
F.tot_ht = F.tot_ht + F.ind_ht;
H.Mc2Source_SPLslider(2) = uicontrol('Style','slider',...
   'Units','pixels',...
   'Position',[0 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd-60 F.ind_ht*.75],...
   'BackgroundColor','red',...
   'Min',0,'Max',30,...
   'SliderStep',[1/31 5/31],...
   'Value',XStimParams.focalSPL(2),...
   'Callback','setinfo_Mc2Source');
H.Mc2Source_SPLslider_txt(2) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[F.Mc2Source_wd-60 F.Mc2Source_ht-F.tot_ht 60 F.ind_ht*.75],...
   'FontSize',8,...
   'String', [num2str(XStimParams.focalSPL(2)) ' dB']);

% SPL inc
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.Mc2Source_ht-F.tot_ht 150 F.ind_ht],...
   'String', 'SPL increment (dB)');
H.Mc2Source_SPLinc = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[155 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
   'String', num2str(XStimParams.SPLinc(1)),...
   'Callback','setinfo_Mc2Source');

% epoch duration
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.Mc2Source_ht-F.tot_ht 150 F.ind_ht],...
   'String', 'epoch DURation (msec)');
H.epoch_duration = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[155 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
   'String', num2str(XStimParams.epoch_duration(1)),...
   'Callback','setinfo_Mc2Source');

% some helpful(?) text
H.Mc2Src_txt(6) = uicontrol('Style','text',...
   'Units','pixels',...
   'fontsize',12,...
   'BackGroundColor', [.8 .8 .8],...
   'Position',[220 F.Mc2Source_ht-F.tot_ht-10 170 F.ind_ht],...
   'String', 'src2 non-focal SPLs');

  % FOCAL total probability
F.tot_ht = F.tot_ht + F.ind_ht;
H.Mc2Src_txt(16) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.Mc2Source_ht-F.tot_ht 150 F.ind_ht],...
   'String', 'source2 FOCAL total PROB');
H.focalProb = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[155 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
   'String', num2str(XStimParams.focalProb),...
   'Callback','setinfo_Mc2Source');

% some helpful(?) text
H.Mc2Src_txt(7) = uicontrol('Style','text',...
   'Units','pixels',...
   'fontsize',12,...
   'BackGroundColor', [.8 .8 .8],...
   'Position',[220 F.Mc2Source_ht-F.tot_ht-10 170 F.ind_ht],...
   'String', 'range from 0 to 30 dB');

    % Source1 SPL
F.tot_ht = F.tot_ht + F.ind_ht;
H.Mc2Src_txt(17) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.Mc2Source_ht-F.tot_ht 150 F.ind_ht],...
   'String', 'source1 SPL (0 to 30dB)');
H.source1_SPL = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[155 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
   'String', num2str(XStimParams.source1_SPL),...
   'Callback','setinfo_Mc2Source');

% some helpful(?) text
H.Mc2Src_txt(8) = uicontrol('Style','text',...
   'Units','pixels',...
   'fontsize',12,...
   'BackGroundColor', [.8 .8 .8],...
   'Position',[210 F.Mc2Source_ht-F.tot_ht-10 180 F.ind_ht],...
   'String', 'with p = 1 - p(Focal)');

% Source1 ISI
F.tot_ht = F.tot_ht + F.ind_ht;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.Mc2Source_ht-F.tot_ht 150 F.ind_ht],...
   'FontWeight','bold',...
   'String', 'source1 ISI (epochs)');
H.source1_ISI = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[155 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
   'String', num2str(XStimParams.source1_ISI),...
   'Callback','setinfo_Mc2Source');
  
% ramp probe?
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[215 F.Mc2Source_ht-F.tot_ht-F.ind_ht*.5 60 F.ind_ht],...
   'String', 'probe ramp?');
H.source1_ramp = uicontrol('Style','checkbox',...
   'Units','pixels',...
   'Position',[275 F.Mc2Source_ht-F.tot_ht-F.ind_ht*.5 15 F.ind_ht],...
   'String', num2str(XStimParams.source1_ISI),...
   'Callback','setinfo_Mc2Source');

% plot stimulus in separate window
H.Mc2Src_plotstim = uicontrol('Style','pushbutton',...
    'Units','pixels',...
    'Position',[320 F.Mc2Source_ht-F.tot_ht-F.ind_ht*.5 80 F.ind_ht],...
    'BackgroundColor','green',...
    'ForegroundColor','blue',...
    'String', 'plot stim',...
    'Callback','Mc2src_plotstim');


%%%% repeat mode controls
F.tot_ht = F.tot_ht + F.ind_ht*1.5;
H.Mc2Src_txt(9) = uicontrol('Style','text',...
    'Units','pixels',...
   'BackGroundColor', [.8 .8 .8],...
   'fontsize',12,...
   'Position',[170 F.Mc2Source_ht-F.tot_ht 195 F.ind_ht],...
    'String', 'total stimdur --> # repeats');
H.Mc2Src_repeatmodepb = uicontrol('Style','togglebutton',...
    'Units','pixels',...
    'Position',[0 F.Mc2Source_ht-F.tot_ht 150 F.ind_ht],...
    'BackgroundColor','yellow',...
    'ForegroundColor','blue',...
    'String', 'push for repeat mode',...
    'Callback','setinfo_Mc2Source');

% on duration    
F.tot_ht = F.tot_ht + F.ind_ht;
H.Mc2Src_txt(10) = uicontrol('Style','text',...
    'Units','pixels',...
    'Position',[0 F.Mc2Source_ht-F.tot_ht 150 F.ind_ht],...
    'String', 'on DURation (#epochs)');
H.Mc2Src_onEpochs = uicontrol('Style','edit',...
    'Units','pixels',...
    'Position',[155 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
    'String', num2str(XStimParams.onEpochs),...
    'Callback','setinfo_Mc2Source');

% some helpful(?) text
H.Mc2Src_txt(11) = uicontrol('Style','text',...
   'Units','pixels',...
   'fontsize',12,...
   'BackGroundColor', [.8 .8 .8],...
   'Position',[205 F.Mc2Source_ht-F.tot_ht 95 F.ind_ht],...
   'String', '(on = n * ISI)');
% probe when adapter on?
H.Mc2Src_txt(14) = uicontrol('Style','text',...
    'Units','pixels',...
    'Position',[300 F.Mc2Source_ht-F.tot_ht 80 F.ind_ht],...
    'BackGroundColor', [.8 .8 .8],...
    'ForegroundColor','blue',...
    'String', 'probe when on?');
H.McSrc_probeFlag(1) = uicontrol('Style','checkbox',...
    'Units','pixels',...
    'Position',[380 F.Mc2Source_ht-F.tot_ht 20 F.ind_ht],...
    'BackGroundColor', [.8 .8 .8],...
    'Value',XStimParams.probeFlag(1),...
    'Callback', 'setinfo_Mc2Source');

% off duration
F.tot_ht = F.tot_ht + F.ind_ht;
H.Mc2Src_txt(12) = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 F.Mc2Source_ht-F.tot_ht 150 F.ind_ht],...
   'String', 'off DURation (#epochs)');
H.Mc2Src_offEpochs = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[155 F.Mc2Source_ht-F.tot_ht 50 F.ind_ht],...
   'String', num2str(XStimParams.offEpochs),...
   'Callback','setinfo_Mc2Source');

% some helpful(?) text
H.Mc2Src_txt(13) = uicontrol('Style','text',...
   'Units','pixels',...
   'BackGroundColor', [.8 .8 .8],...
   'fontsize',12,...
   'Position',[205 F.Mc2Source_ht-F.tot_ht 95 F.ind_ht],...
   'String', '(off = n * ISI)');
    
% probe when adapter off?
H.Mc2Src_txt(15) = uicontrol('Style','text',...
    'Units','pixels',...
    'Position',[300 F.Mc2Source_ht-F.tot_ht 80 F.ind_ht],...
    'BackGroundColor', [.8 .8 .8],...
    'ForegroundColor','blue',...
    'String', 'probe when off?');
H.McSrc_probeFlag(2) = uicontrol('Style','checkbox',...
    'Units','pixels',...
    'Position',[380 F.Mc2Source_ht-F.tot_ht 20 F.ind_ht],...
    'BackGroundColor', [.8 .8 .8],...
    'value',XStimParams.probeFlag(2),...
    'Callback', 'setinfo_Mc2Source');
    
    
    
    % run control
    F.tot_ht = F.tot_ht + F.ind_ht*1.5;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Run Control');
    
    %Engage TEST
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.engageMc2Source = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','green',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'ENGAGE',...
        'Callback','Engage_Mc2Source_new');
    
    %Record Data?
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[85 F.Mc2Source_ht-F.tot_ht 100 F.ind_ht],...
        'String', 'Record Data?');
    H.recorddata = uicontrol('Style','checkbox',...
        'Units','pixels',...
        'Position',[185 F.Mc2Source_ht-F.tot_ht 20 F.ind_ht],...
        'Callback', 'setinfo_Mc2Source');
    F.tot_ht = F.tot_ht + F.ind_ht;
    H.recorddata_FN = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[85 F.Mc2Source_ht-F.tot_ht 100 F.ind_ht],...
        'String', FN.data);
    
    % increment test number
    H.inc_testnum = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[250 F.Mc2Source_ht-F.tot_ht 140 F.ind_ht],...
        'BackgroundColor','yellow',...
        'ForegroundColor','blue',...
        'FontWeight','bold',...
        'String', 'INCREMENT testnum',...
        'Callback','setInfo_Mc2Source');
    
    %Pause TEST
    F.ind_ht = 20; 
    F.tot_ht = F.tot_ht + F.ind_ht + F.pos_incr;
    H.pauseMc2Source = uicontrol('Style','togglebutton',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Pause');
   
    %Reset
    H.resetMc2Source = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[200 F.Mc2Source_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','cyan',...
        'ForegroundColor','black',...
        'FontWeight','bold',...
        'String', 'Reset',...
        'Callback','Reset_Mc2Source');
    
    %Exit
    H.exitMc2Source = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[F.Mc2Source_wd-100 F.Mc2Source_ht-F.tot_ht 80 F.ind_ht],...
        'BackgroundColor','red',...
        'ForegroundColor','yellow',...
        'FontWeight','bold',...
        'String', 'EXIT',...
        'Callback','Exit_Mc2Source');
    
    %Build/Play Trials
    F.ind_ht = 20; 
    F.tot_ht = F.tot_ht + F.ind_ht + F.pos_incr;
    H.status = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht F.Mc2Source_wd F.ind_ht],...
        'BackgroundColor','blue',...
        'ForegroundColor','white',...
        'FontWeight','bold',...
        'String', 'Status: setting parameters');
    
    %Remaining Reps and Trials
    F.tot_ht = F.tot_ht + F.ind_ht;
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[0 F.Mc2Source_ht-F.tot_ht 120 F.ind_ht],...
        'String', 'Remaining reps:');
    H.remreps = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[122 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
        'String', get(H.numreps,'String'));
    uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[164 F.Mc2Source_ht-F.tot_ht 80 F.ind_ht],...
        'String', 'Trials:');
    H.remtrials = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[246 F.Mc2Source_ht-F.tot_ht 40 F.ind_ht],...
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

SetInfo_Mc2Source;