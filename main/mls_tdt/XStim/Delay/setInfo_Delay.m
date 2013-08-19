% SetInfo_Delay
% combines all the other 'set...' scripts into one


global XStimParams;

XStimParams.space_flag = get(H.spaceflag,'Value');
%disp(XStimParams.space_flag);
%set(H.spaceflag,'Value',XStimParams.space_flag);

%control visibility of file picking for ILA, ITA or fully-cued
%control visibility of file picking for ILA, ITA or fully-cued
if XStimParams.space_flag
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');   
else
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');   
end


% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
    [FN.space,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1) = testHRTFfiletype(FN.space_path, FN.space);
end

% set reference locations
% if get(H.spaceflag,'Value')
%     XStimParams.az = str2num(get(H.az,'String'));
%     XStimParams.el = str2num(get(H.el,'String'));
% end

% set numreps
XStimParams.numreps = str2num(get(H.numreps,'String'));

% set offsets
XStimParams.offset_el(1) = 0;
XStimParams.offset_az(1) = 0;
XStimParams.offset_el(2) = str2num(get(H.offset_el,'String'));
XStimParams.offset_az(2) = str2num(get(H.offset_az,'String'));


% set multiplication factor for each source
% XStimParams.factor(1) = str2num(get(H.Delay_factor,'String'));
% XStimParams.factor(2) = str2num(get(H.Delay_factor2,'String'));


%Set space PickerFlag
XStimParams.picker_flag = 1;
if exist1('H.pickerfig') 	
    if isempty(H.pickerfig)
        XStimParams.locations = [];
        H.pickerfig = spacePicker;
    end
else
    XStimParams.locations = [];
    H.pickerfig = spacePicker;  
end

% increment test number
if(exist1('H.Delayfig') & get(H.inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.inc_testnum,'Value',0);
end

% update stim_type1
%XStimParams.stim_type = get(H.Delay_stim_type,'String');
%H.stim_type = get(H.stim_type,'Value');

%control visibility of FN vs freq
%set(H.Delay_stim_filenamepb,'Visible','on');
% set(H.Delay_stim_filename,'Visible','on');
% set(H.Delay_freq_txt,'Visible','off');
% set(H.Delay_freq,'Visible','off');    

% switch deblank(XStimParams.stim_type)
%     case 'File'
%         set(H.Delay_stim_filename,'Enable','on');
%         set(H.stim_filename,'Enable','on');
%     case {'Gammatone', 'Tone'}
%         set(H.Delay_freq_txt,'Visible','on');
%         set(H.Delay_freq,'Visible','on');    
%         %set(H.Delay_stim_filenamepb,'Visible','off');
%         set(H.Delay_stim_filename,'Visible','off');
%     otherwise
%         set(H.Delay_stim_filename,'Enable','off');
%         set(H.stim_filename,'Enable','off');
% end


% set dur based on stim1
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));


% % bandwidths and bandlimits
% % 1st stim
% if ~strcmp(XStimParams.stim_type, 'File') 	% if not from file
%     % if changed, re-set bandwidths
%     if XStimParams.bandwidth(1) ~= str2num(get(H.bandwidth1,'String'))
%         XStimParams.bandwidth(1) = str2num(get(H.bandwidth1,'String'));
%     end
%     % if changed, re-set bandlimits
%     if (XStimParams.bandlimits(1) ~= str2num(get(H.bandlimits1,'String')) | ...
%             XStimParams.bandlimits(2) ~= str2num(get(H.bandlimits2,'String')));  
%         XStimParams.bandlimits(1) = str2num(get(H.bandlimits1,'String'));
%         XStimParams.bandlimits(2) = str2num(get(H.bandlimits2,'String'));
%     end
%     set(H.bandwidth1,'Enable','on');
%     set(H.bandlimits1,'Enable','on');
%     set(H.bandlimits2,'Enable','on');
% else
%     set(H.bandwidth1,'Enable','off');
%     set(H.bandlimits1,'Enable','off');
%     set(H.bandlimits2,'Enable','off');
% end


XStimParams.test_ISI = str2num(get(H.ISI,'string'));

% Delay and Onset Offset variables
XStimParams.randOnsetPerms = get(H.randOnsetPerms,'Value');

% - XStimParams.Delay_ms = str2num(get(H.delaytime,'string'));
XStimParams.uncorrel = get(H.uncorrel,'Value');
if(XStimParams.uncorrel & XStimParams.SilentDelay & ~XStimParams.randOnsetPerms)
    XStimParams.SilentDelay = 0;
    set(H.silentdelay, 'Value', XStimParams.SilentDelay);
end

XStimParams.SilentDelay = get(H.silentdelay,'Value');
if(XStimParams.SilentLag & XStimParams.SilentDelay & ~XStimParams.randOnsetPerms)
    XStimParams.SilentLag = 0;
    set(H.silentlag, 'Value', XStimParams.SilentLag);
end
if(XStimParams.uncorrel & XStimParams.SilentDelay & ~XStimParams.randOnsetPerms)
    XStimParams.uncorrel = 0;
    set(H.uncorrel, 'Value', XStimParams.uncorrel);
end

XStimParams.SilentLag = get(H.silentlag,'Value');

if(XStimParams.SilentLag & XStimParams.SilentDelay & ~XStimParams.randOnsetPerms)
    XStimParams.SilentDelay = 0;
    set(H.silentdelay, 'Value', XStimParams.SilentDelay);
end
if(XStimParams.uncorrel & XStimParams.SilentLag & ~XStimParams.randOnsetPerms)
    XStimParams.uncorrel = 0;
    set(H.uncorrel, 'Value', XStimParams.uncorrel);
end


% onset offset ramps
% - XStimParams.ramp_time_on = str2num(get(H.onramptime,'string'));
% - XStimParams.ramp_time_off = str2num(get(H.offramptime,'string'));

XStimParams.DelayTimes(1) = str2num(get(H.delay01,'String'));
XStimParams.DelayTimes(2) = str2num(get(H.delay02,'String'));
XStimParams.DelayTimes(3) = str2num(get(H.delay03,'String'));
XStimParams.DelayTimes(4) = str2num(get(H.delay04,'String'));
XStimParams.DelayTimes(5) = str2num(get(H.delay05,'String'));
XStimParams.DelayTimes(6) = str2num(get(H.delay06,'String'));
XStimParams.DelayTimes(7) = str2num(get(H.delay07,'String'));
XStimParams.DelayTimes(8) = str2num(get(H.delay08,'String'));
XStimParams.DelayTimes(9) = str2num(get(H.delay09,'String'));
XStimParams.DelayTimes(10) = str2num(get(H.delay10,'String'));
XStimParams.DelayTimes(11) = str2num(get(H.delay11,'String'));
XStimParams.DelayTimes(12) = str2num(get(H.delay12,'String'));
XStimParams.DelayTimes(13) = str2num(get(H.delay13,'String'));
XStimParams.DelayTimes(14) = str2num(get(H.delay14,'String'));
XStimParams.DelayTimes(15) = str2num(get(H.delay15,'String'));
XStimParams.DelayTimes(16) = str2num(get(H.delay16,'String'));

XStimParams.ramp_timeS_on(1) = str2num(get(H.onramp01,'String'));
XStimParams.ramp_timeS_on(2) = str2num(get(H.onramp02,'String'));
XStimParams.ramp_timeS_on(3) = str2num(get(H.onramp03,'String'));
XStimParams.ramp_timeS_on(4) = str2num(get(H.onramp04,'String'));
XStimParams.ramp_timeS_on(5) = str2num(get(H.onramp05,'String'));
XStimParams.ramp_timeS_on(6) = str2num(get(H.onramp06,'String'));

XStimParams.ramp_timeS_off(1) = str2num(get(H.offramp01,'String'));
XStimParams.ramp_timeS_off(2) = str2num(get(H.offramp02,'String'));
XStimParams.ramp_timeS_off(3) = str2num(get(H.offramp03,'String'));
XStimParams.ramp_timeS_off(4) = str2num(get(H.offramp04,'String'));
XStimParams.ramp_timeS_off(5) = str2num(get(H.offramp05,'String'));
XStimParams.ramp_timeS_off(6) = str2num(get(H.offramp06,'String'));

XStimParams.lag_seg_ext(1) = str2num(get(H.lagext01,'String'));
XStimParams.lag_seg_ext(2) = str2num(get(H.lagext02,'String'));
XStimParams.lag_seg_ext(3) = str2num(get(H.lagext03,'String'));
XStimParams.lag_seg_ext(4) = str2num(get(H.lagext04,'String'));
XStimParams.lag_seg_ext(5) = str2num(get(H.lagext05,'String'));
XStimParams.lag_seg_ext(6) = str2num(get(H.lagext06,'String'));
XStimParams.lag_seg_ext(7) = str2num(get(H.lagext07,'String'));
XStimParams.lag_seg_ext(8) = str2num(get(H.lagext08,'String'));



if(XStimParams.randOnsetPerms)
    set(H.delay01,		    'Enable', 'on');
    set(H.delay02,		    'Enable', 'on');
    set(H.delay03,		    'Enable', 'on');
    set(H.delay04,		    'Enable', 'on');
    set(H.delay05,		    'Enable', 'on');
    set(H.delay06,		    'Enable', 'on');
    set(H.delay07,		    'Enable', 'on');
    set(H.delay08,		    'Enable', 'on');
    set(H.delay09,		    'Enable', 'on');
    set(H.delay10,		    'Enable', 'on');
    set(H.delay11,		    'Enable', 'on');
    set(H.delay12,		    'Enable', 'on');
    set(H.delay13,		    'Enable', 'on');
    set(H.delay14,		    'Enable', 'on');
    set(H.delay15,		    'Enable', 'on');
    set(H.delay16,		    'Enable', 'on');

    set(H.DelayOnOff_mode,		'Enable', 'on');
    
    %set(H.onramp01,		    'Enable', 'on');
    set(H.onramp02,		    'Enable', 'on');
    set(H.onramp03,		    'Enable', 'on');
    set(H.onramp04,		    'Enable', 'on');
    set(H.onramp05,		    'Enable', 'on');
    set(H.onramp06,		    'Enable', 'on');
    
    %set(H.offramp01,		    'Enable', 'on');
    set(H.offramp02,		    'Enable', 'on');
    set(H.offramp03,		    'Enable', 'on');
    set(H.offramp04,		    'Enable', 'on');
    set(H.offramp05,		    'Enable', 'on');
    set(H.offramp06,		    'Enable', 'on');

    % make sure onsets and odffsets are positive (window not updated)
    for d=1:size(XStimParams.ramp_timeS_on,2)
        if(XStimParams.ramp_timeS_on(d) < 0)
            XStimParams.ramp_timeS_on(d) = -XStimParams.ramp_timeS_on(d);
            eval(['set(H.onramp0' round(num2str(d)) ', ''String'', ' '''' num2str(XStimParams.ramp_timeS_on(d)) ''');']);
        end
    end
    for d=1:size(XStimParams.ramp_timeS_off,2)
        if(XStimParams.ramp_timeS_off(d) < 0)
            XStimParams.ramp_timeS_off(d) = -XStimParams.ramp_timeS_off(d);
            eval(['set(H.offramp0' round(num2str(d)) ', ''String'', ' '''' num2str(XStimParams.ramp_timeS_off(d)) ''');']);
        end
    end
    
    set(H.lagext01,		    'Enable', 'on');
    set(H.lagext02,		    'Enable', 'on');
    set(H.lagext03,		    'Enable', 'on');
    set(H.lagext04,		    'Enable', 'on');
    set(H.lagext05,		    'Enable', 'on');
    set(H.lagext06,		    'Enable', 'on');
    set(H.lagext07,		    'Enable', 'on');
    set(H.lagext08,		    'Enable', 'on');
    
    if(XStimParams.lag_seg_ext(1)>0)
        warning('The first lag extension should probably be 0');    
    end
    
else
    % - set(H.delaytime,		'Enable', 'on');
    set(H.uncorrel,		    'Enable', 'on');
    set(H.silentdelay,		'Enable', 'on');
    
    % - set(H.onramptime,		'Enable', 'on');
    % - set(H.offramptime,		'Enable', 'on');
    set(H.DelayOnOff_mode,  'Enable', 'off');
    
    %set(H.delay01,		    'Enable', 'off');
    set(H.delay02,		    'Enable', 'off');
    set(H.delay03,		    'Enable', 'off');
    set(H.delay04,		    'Enable', 'off');
    set(H.delay05,		    'Enable', 'off');
    set(H.delay06,		    'Enable', 'off');
    set(H.delay07,		    'Enable', 'off');
    set(H.delay08,		    'Enable', 'off');
    set(H.delay09,		    'Enable', 'off');
    set(H.delay10,		    'Enable', 'off');
    set(H.delay11,		    'Enable', 'off');
    set(H.delay12,		    'Enable', 'off');
    set(H.delay13,		    'Enable', 'off');
    set(H.delay14,		    'Enable', 'off');
    set(H.delay15,		    'Enable', 'off');
    set(H.delay16,		    'Enable', 'off');
    
    %set(H.onramp01,		    'Enable', 'off');
    set(H.onramp02,		    'Enable', 'off');
    set(H.onramp03,		    'Enable', 'off');
    set(H.onramp04,		    'Enable', 'off');
    set(H.onramp05,		    'Enable', 'off');
    set(H.onramp06,		    'Enable', 'off');
    %set(H.onramp07,		    'Enable', 'off');
    %set(H.onramp08,		    'Enable', 'off');
   % set(H.onramp09,		    'Enable', 'off');
    
    %set(H.offramp01,		    'Enable', 'off');
    set(H.offramp02,		    'Enable', 'off');
    set(H.offramp03,		    'Enable', 'off');
    set(H.offramp04,		    'Enable', 'off');
    set(H.offramp05,		    'Enable', 'off');
    set(H.offramp06,		    'Enable', 'off');
    %set(H.offramp07,		    'Enable', 'off');
    %set(H.offramp08,		    'Enable', 'off');
    %set(H.offramp09,		    'Enable', 'off');
    
    %set(H.lagext01,		    'Enable', 'off');
    set(H.lagext02,		    'Enable', 'off');
    set(H.lagext03,		    'Enable', 'off');
    set(H.lagext04,		    'Enable', 'off');
    set(H.lagext05,		    'Enable', 'off');
    set(H.lagext06,		    'Enable', 'off');
    set(H.lagext07,		    'Enable', 'off');
    set(H.lagext08,		    'Enable', 'off');

end

XStimParams.curr_ABL = str2num(get(H.ABL,'string'));

%check visibility of record datafile
if get(H.recorddata,'Value')
    set(H.recorddata_FN,'enable','on');
else
    set(H.recorddata_FN,'Enable','off');
end

% % check out modulation parameters for stim1
% XStimParams.mod_depth(1) = str2num(get(H.Delay_mod_depth,'String'));
% mod_type = get(H.Delay_mod_type,'String');
% mod_num = get(H.Delay_mod_type,'Value');
% XStimParams.mod_type = deblank(mod_type(mod_num,:));


% switch XStimParams.mod_type
%     case 'Tone'
%         set(H.Delay_mod_txt,'String', 'Mod Freq:');
%         set(H.Delay_mod_freq,'Visible','on');
%         set(H.Delay_mod_phase,'Visible','on');
%         set(H.Delay_mod_txtA,'Visible','on');
%         set(H.Delay_mod_depth,'Visible','on');
%         set(H.Delay_mod_txtB,'Visible','on');
%         set(H.Delay_mod_pb,'Visible','off');
%         XStimParams.mod_freq(1) = str2num(get(H.Delay_mod_freq,'String'));
%         XStimParams.mod_phase(1) = str2num(get(H.Delay_mod_phase,'String'));
%     case 'File'
%         set(H.Delay_mod_pb,'Visible','on');
%         set(H.Delay_mod_freq,'Visible','off');
%         set(H.Delay_mod_txtA,'Visible','on');
%         set(H.Delay_mod_depth,'Visible','on');
%         set(H.Delay_mod_phase,'Visible','off');
%         set(H.Delay_mod_txtB,'Visible','off');
%         if get(H.Delay_mod_pb,'Value')
%             [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
%             set(H.Delay_mod_pb,'Value',0);
%         end
%         set(H.Delay_mod_txt,'String', FN.mod);
%         
%     case 'LP Noise'
%         set(H.Delay_mod_txt,'String', 'CutOff Freq:');
%         set(H.Delay_mod_freq,'Visible','on');
%         set(H.Delay_mod_txtA,'Visible','on');
%         set(H.Delay_mod_depth,'Visible','on');
%         set(H.Delay_mod_phase,'Visible','off');
%         set(H.Delay_mod_pb,'Visible','off');
%         set(H.Delay_mod_txtB,'Visible','off');
%         XStimParams.mod_freq(1) = str2num(get(H.Delay_mod_freq,'String'));
%         
%     case 'None'
%         set(H.Delay_mod_txt,'String', 'no mod  ');
%         set(H.Delay_mod_freq,'Visible','off');
%         set(H.Delay_mod_phase,'Visible','off');
%         set(H.Delay_mod_pb,'Visible','off');
%         set(H.Delay_mod_txtB,'Visible','off');
%         set(H.Delay_mod_txtA,'Visible','off');
%         set(H.Delay_mod_depth,'Visible','off');
%         
%     otherwise
% end

% % check out modulation parameters for stim2
% XStimParams.mod_depth(2) = str2num(get(H.Delay_mod_depth2,'String'));
% mod_type = get(H.Delay_mod_type2,'String');
% mod_num = get(H.Delay_mod_type2,'Value');
% XStimParams.mod_type2 = deblank(mod_type(mod_num,:));
% 
% switch XStimParams.mod_type2
%     case 'Tone'
%         set(H.Delay_mod_txt2,'String', 'Mod Freq:');
%         set(H.Delay_mod_freq2,'Visible','on');
%         set(H.Delay_mod_pb2,'Visible','off');
%         set(H.Delay_mod_txtA2,'Visible','on');
%         set(H.Delay_mod_depth2,'Visible','on');
%         set(H.Delay_mod_phase2,'Visible','on');
%         set(H.Delay_mod_txtB2,'Visible','on');
%         XStimParams.mod_freq(2) = str2num(get(H.Delay_mod_freq2,'String'));
%         XStimParams.mod_phase(2) = str2num(get(H.Delay_mod_phase2,'String'));
%     case 'File'
%         set(H.Delay_mod_pb2,'Visible','on');
%         set(H.Delay_mod_freq2,'Visible','off');
%         set(H.Delay_mod_txtA2,'Visible','on');
%         set(H.Delay_mod_depth2,'Visible','on');
%         set(H.Delay_mod_phase2,'Visible','off');
%         set(H.Delay_mod_txtB2,'Visible','off');
%         if get(H.Delay_mod_pb2,'Value')
%             [mod_from_file2, FN.mod2, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
%             set(H.Delay_mod_pb2,'Value',0);
%         end
%         set(H.Delay_mod_txt2,'String', FN.mod2);
%         
%     case 'LP Noise'
%         set(H.Delay_mod_txt2,'String', 'CutOff Freq:');
%         set(H.Delay_mod_freq2,'Visible','on');
%         set(H.Delay_mod_txtA2,'Visible','on');
%         set(H.Delay_mod_depth2,'Visible','on');
%         set(H.Delay_mod_pb2,'Visible','off');
%         set(H.Delay_mod_phase2,'Visible','off');
%         set(H.Delay_mod_txtB2,'Visible','off');
%         XStimParams.mod_freq(2) = str2num(get(H.Delay_mod_freq2,'String'));
%         
%     case 'None'
%         set(H.Delay_mod_txt2,'String', 'no mod  ');
%         set(H.Delay_mod_freq2,'Visible','off');
%         set(H.Delay_mod_phase2,'Visible','off');
%         set(H.Delay_mod_pb2,'Visible','off');
%         set(H.Delay_mod_txtB2,'Visible','off');
%         set(H.Delay_mod_txtA2,'Visible','off');
%         set(H.Delay_mod_depth2,'Visible','off');
%         
%     otherwise
% end
% 
% % check out modulation parameters for stim3
% XStimParams.mod_depth(3) = str2num(get(H.Delay_mod_depth3,'String'));
% mod_type = get(H.Delay_mod_type3,'String');
% mod_num = get(H.Delay_mod_type3,'Value');
% XStimParams.mod_type3 = deblank(mod_type(mod_num,:));
% 
% switch XStimParams.mod_type3
%     case 'Tone'
%         set(H.Delay_mod_txt3,'String', 'Mod Freq:');
%         set(H.Delay_mod_freq3,'Visible','on');
%         set(H.Delay_mod_txtA3,'Visible','on');
%         set(H.Delay_mod_depth3,'Visible','on');
%         set(H.Delay_mod_phase3,'Visible','on');
%         set(H.Delay_mod_pb3,'Visible','off');
%         set(H.Delay_mod_txtB3,'Visible','on');
%         XStimParams.mod_freq(3) = str2num(get(H.Delay_mod_freq3,'String'));
%         XStimParams.mod_phase(3) = str2num(get(H.Delay_mod_phase3,'String'));
%     case 'File'
%         set(H.Delay_mod_pb3,'Visible','on');
%         set(H.Delay_mod_freq3,'Visible','off');
%         set(H.Delay_mod_txtA3,'Visible','on');
%         set(H.Delay_mod_depth3,'Visible','on');
%         set(H.Delay_mod_phase3,'Visible','off');
%         set(H.Delay_mod_txtB3,'Visible','off');
%         if get(H.Delay_mod_pb3,'Value')
%             [mod_from_file3, FN.mod3, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
%             set(H.Delay_mod_pb3,'Value',0);
%         end
%         set(H.Delay_mod_txt3,'String', FN.mod2);
%         
%     case 'LP Noise'
%         set(H.Delay_mod_txt3,'String', 'CutOff Freq:');
%         set(H.Delay_mod_freq3,'Visible','on');
%         set(H.Delay_mod_txtA3,'Visible','on');
%         set(H.Delay_mod_depth3,'Visible','on');
%         set(H.Delay_mod_pb3,'Visible','off');
%         set(H.Delay_mod_phase3,'Visible','off');
%         set(H.Delay_mod_txtB3,'Visible','off');
%         XStimParams.mod_freq(3) = str2num(get(H.Delay_mod_freq3,'String'));
%         
%     case 'None'
%         set(H.Delay_mod_txt3,'String', 'no mod  ');
%         set(H.Delay_mod_freq3,'Visible','off');
%         set(H.Delay_mod_pb3,'Visible','off');
%         set(H.Delay_mod_phase3,'Visible','off');
%         set(H.Delay_mod_txtB3,'Visible','off');
%         set(H.Delay_mod_txtA3,'Visible','off');
%         set(H.Delay_mod_depth3,'Visible','off');
%         
%     otherwise
% end


% % calc Vstrength
% if get(H.plotVstr,'value')
%     tit = 'Vector Strength Plots';
%     prompt = {'Enter datapath', 'Enter filename ','Enter vector of modulation frequencies'};
%     DefAns = {FN.data_path,[FN.data '.mat'], '[55 75]'};
%     
%     Ans = inputdlg(prompt,tit,1,DefAns);
%     read_data3_minus_spont(Ans{2},Ans{1}, [], 1,str2num(Ans{3}));
%     clear tit prompt DefAns Ans
%     set(H.plotVstr,'value',0);
% end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
