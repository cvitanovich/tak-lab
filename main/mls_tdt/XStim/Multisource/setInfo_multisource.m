% setInfo_multisource
% combines all the other 'set...' scripts into one

if get(H.spaceflag,'Value') & (XStimParams.space_flag == 0)
    XStimParams.space_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
    XStimParams.ABLalone_flag = 0;
elseif get(H.ABLaloneflag,'Value') & (XStimParams.ABLalone_flag == 0)
    XStimParams.space_flag = 0;
    XStimParams.ABLalone_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
elseif get(H.ildaloneflag,'Value') & (XStimParams.ildalone_flag == 0)
    XStimParams.space_flag = 0;
    XStimParams.ildalone_flag = 1;
    XStimParams.itdalone_flag = 0;
    XStimParams.ABLalone_flag = 0;
elseif get(H.itdaloneflag,'Value') & (XStimParams.itdalone_flag == 0)
    XStimParams.space_flag = 0;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 1;
    XStimParams.ABLalone_flag = 0;
end
set(H.spaceflag,'Value',XStimParams.space_flag);
set(H.ildaloneflag,'value',XStimParams.ildalone_flag);
set(H.itdaloneflag,'value',XStimParams.itdalone_flag);
set(H.ABLaloneflag,'value',XStimParams.ABLalone_flag);

% control visibility of file picking for ILA, ITA or fully-cued
if XStimParams.space_flag
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');   
    set(H.ildalonefilehdr,'Enable','off');
    set(H.ildalonefile,'Enable','off');
    set(H.ildalonepb,'Enable','off');  
    set(H.itdalonefilehdr,'Enable','off');
    set(H.itdalonefile,'Enable','off');
    set(H.itdalonepb,'Enable','off');  
    set(H.ITDtxt,'BackgroundColor',[.8 .8 .8]);
    set(H.ILDtxt,'BackgroundColor',[.8 .8 .8]);
elseif XStimParams.ABLalone_flag
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');   
    set(H.ildalonefilehdr,'Enable','off');
    set(H.ildalonefile,'Enable','off');
    set(H.ildalonepb,'Enable','off');  
    set(H.itdalonefilehdr,'Enable','off');
    set(H.itdalonefile,'Enable','off');
    set(H.itdalonepb,'Enable','off');  
    set(H.ITDtxt,'BackgroundColor','green');
    set(H.ILDtxt,'BackgroundColor','green');
elseif XStimParams.ildalone_flag
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');
    set(H.ildalonefilehdr,'Enable','on');
    set(H.ildalonefile,'Enable','on');
    set(H.ildalonepb,'Enable','on');
    set(H.itdalonefilehdr,'Enable','off');
    set(H.itdalonefile,'Enable','off');
    set(H.itdalonepb,'Enable','off');  
    set(H.ITDtxt,'BackgroundColor','green');
    set(H.ILDtxt,'BackgroundColor',[.8 .8 .8]);
elseif XStimParams.itdalone_flag
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');
    set(H.ildalonefilehdr,'Enable','off');
    set(H.ildalonefile,'Enable','off');
    set(H.ildalonepb,'Enable','off');
    set(H.itdalonefilehdr,'Enable','on');
    set(H.itdalonefile,'Enable','on');
    set(H.itdalonepb,'Enable','on');  
    set(H.ITDtxt,'BackgroundColor',[.8 .8 .8]);
    set(H.ILDtxt,'BackgroundColor','green');
end

% check for IA pushbutton (file-picker)
if get(H.ildalonepb,'Value')
    [FN.ildalone,FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select ILD-alone HRTF File');
    if(FN.ILA_path ~= 0)
        set(H.ildalonefile,'String',FN.ildalone);
    end
    set(H.ildalonepb,'Value',0);
    FN.HRTFfiletype(2) = testHRTFfiletype(FN.ILA_path, FN.ildalone);
end

% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
    [FN.space,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select Fully-cued HRTF File');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1) = testHRTFfiletype(FN.space_path, FN.space);
end

% check for ITA pushbutton (file-picker)
if get(H.itdalonepb,'Value')
    [FN.itdalone,FN.ITA_path] = uigetfile([FN.ITA_path '*.*'],'Select ITD-alone HRTF File');
    if(FN.ITA_path ~= 0)
        set(H.itdalonefile,'String',FN.itdalone);
    end
    set(H.itdalonepb,'Value',0);
    FN.HRTFfiletype(3) = testHRTFfiletype(FN.ITA_path, FN.itdalone);
end

% set reference locations
if get(H.spaceflag,'Value') | get(H.ABLaloneflag,'Value') | get(H.itdaloneflag,'Value')
    XStimParams.az = str2num(get(H.az,'String'));
    XStimParams.el = str2num(get(H.el,'String'));
end

% set numreps
XStimParams.numreps = str2num(get(H.numreps,'String'));

% increment test number
if(exist1('H.multi_sourcefig') & get(H.inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.inc_testnum,'Value',0);
end

% % % update stim_type
% % stim_type = get(H.Twosrc_stim_type,'String');
% % stim_val = get(H.Twosrc_stim_type,'Value');
% % XStimParams.stim_type = deblank(stim_type(stim_val,:));
% % 
% % %if file pb pushed
% % if ~strcmp(XStimParams.stim_type, deblank(stim_type(stim_val,:))) | get(H.Twosrc_stim_filenamepb,'Value')
% %     set(H.stim_type,'Value',stim_val);
% %     if(strcmp('File',deblank(XStimParams.stim_type)))
% %         [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
% %         set(H.Twosrc_stim_filename,'String', FN.stim);
% %         set(H.stim_filename,'String', FN.stim);
% %     end
% %     %setinfo_Xstim;
% %     set(H.Twosrc_stim_filenamepb,'Value',0);
% % end
% 
% %control visibility of FN vs freq
% switch deblank(XStimParams.stim_type)
%     case {'Gammatone', 'Tone'}
%         set(H.Twosrc_stim_filename,'Enable','off');
%         set(H.stim_filename,'Enable','off');
% 		set(H.Twosrc_stim_filenamepb,'Visible','off');
%         set(H.Twosrc_stim_filename,'Visible','off');
% 		set(H.Twosrc_freq_txt,'Visible','on');
%         set(H.Twosrc_freq,'Visible','on');       
% 		set(H.stim_nostim_text,'Visible','off');
% 	case 'File'
%         set(H.Twosrc_stim_filename,'Enable','on');
%         set(H.stim_filename,'Enable','on');
% 		set(H.Twosrc_stim_filenamepb,'Visible','on');
%         set(H.Twosrc_stim_filename,'Visible','on');
% 		set(H.Twosrc_freq_txt,'Visible','off');
%         set(H.Twosrc_freq,'Visible','off');   
% 		set(H.stim_nostim_text,'Visible','off');
% 	case 'None'
% 		set(H.Twosrc_stim_filename,'Enable','off');
%         set(H.stim_filename,'Enable','off');
% 		set(H.Twosrc_stim_filenamepb,'Visible','off');
%         set(H.Twosrc_stim_filename,'Visible','off');
% 		set(H.Twosrc_freq_txt,'Visible','off');
%         set(H.Twosrc_freq,'Visible','off');
% 		set(H.stim_nostim_text,'Visible','on');
% 	otherwise
% 	 	set(H.Twosrc_stim_filename,'Enable','off');
%         set(H.stim_filename,'Enable','off');
% 		set(H.Twosrc_stim_filenamepb,'Visible','off');
%         set(H.Twosrc_stim_filename,'Visible','off');
% 		set(H.Twosrc_freq_txt,'Visible','off');
%         set(H.Twosrc_freq,'Visible','off');
% 		set(H.stim_nostim_text,'Visible','off');
% end

% % buttons to open other windows
% if strcmp(XStimParams.stim_type, 'None')
% 	set(H.source_sub1_gui_pb,'Enable','off');
%     set(H.source_sub2_gui_pb,'Enable','off');
%     set(H.source_sub3_gui_pb,'Enable','off');   
% elseif strcmp(XStimParams.stim_type, 'File') | strcmp(XStimParams.stim_type, 'Gammatone') | strcmp(XStimParams.stim_type, 'Tone')
% 	set(H.source_sub1_gui_pb,'Enable','on');
%     set(H.source_sub2_gui_pb,'Enable','on');
%     set(H.source_sub3_gui_pb,'Enable','off');   
% else
% 	set(H.source_sub1_gui_pb,'Enable','on');
%     set(H.source_sub2_gui_pb,'Enable','on');
%     set(H.source_sub3_gui_pb,'Enable','on');   
% end

if get(H.source_sub1_gui_pb, 'Value')
	H.shape_flag = 1;
	spacePicker_ted;
	set(H.source_sub1_gui_pb, 'Value', 0);
end

if get(H.source_sub2_gui_pb, 'Value')
	H.shape_flag = 2;
	spacePicker_ted;
	set(H.source_sub2_gui_pb, 'Value', 0);
end

if get(H.source_sub3_gui_pb,'Value')
  	clear locations;
	eval(['load ' FN.current_path 'Locations_current1;']);
	XStimParams.sub_locations = locations;
	clear locations;
	for count1 = 1:size(XStimParams.sub_locations, 1)
        eval('XStimParams.stim_type_ted(count1);', 'XStimParams.stim_type_ted(count1) = 8;');
        eval('XStimParams.mod_type_ted(count1);', 'XStimParams.mod_type_ted(count1) = 4;');
        eval('XStimParams.freq(count1);', 'XStimParams.freq(count1) = 100;');
        eval('XStimParams.factor(count1);', 'XStimParams.factor(count1) = 1;');
        eval('XStimParams.bandwidth(count1);', 'XStimParams.bandwidth(count1) = 1000;');
        eval('XStimParams.bandlimits1(count1);', 'XStimParams.bandlimits1(count1) = 10;');
        eval('XStimParams.bandlimits2(count1);', 'XStimParams.bandlimits2(count1) = 1010;');
        eval('XStimParams.mod_freq(count1);', 'XStimParams.mod_freq(count1) = 100;');
        eval('XStimParams.mod_depth(count1);', 'XStimParams.mod_depth(count1) = 1;');
        eval('XStimParams.mod_phase(count1);', 'XStimParams.mod_phase(count1) = 1;');
        H.stim_type(count1) = XStimParams.stim_type_ted(count1);
		H.mod_type(count1) = XStimParams.mod_type_ted(count1);		
        ['stim', num2str(count1), ' el/az ', num2str(XStimParams.sub_locations(count1, :))]
    end
    setParams;
	set(H.source_sub3_gui_pb, 'Value', 0);
end

% if get(H.source_sub3_gui_pb,'Value')
% 	clear locations;
% 	eval(['load ' FN.current_path 'Locations_current1;'])
% 	XStimParams.sub_locations = locations;
% 	clear locations;
% 	for count1 = 1:size(XStimParams.sub_locations, 1)
% 		XStimParams.factor(count1) = 1;
% 	end
% 	setWeights;
% 	set(H.source_sub3_gui_pb,'Value',0);
% end
% 
% if get(H.source_sub4_gui_pb,'Value')
% 	clear locations;
% 	eval(['load ' FN.current_path 'Locations_current1;'])
% 	XStimParams.sub_locations = locations;
% 	clear locations;
% 	for count1 = 1:size(XStimParams.sub_locations, 1)
% 		XStimParams.bandwidth(count1) = 100;
% 		XStimParams.bandlimits1(count1) = 1000;
% 		XStimParams.bandlimits2(count1) = 10000;
% 	end
% 	setNoise;
% 	set(H.source_sub4_gui_pb,'Value',0);
% end
% 
% 
% % set dur 
% XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
% 
% isthereafile_flag = 0;
% for count1 = 1:size(XStimParams.sub_locations, 1)
%     if exist(H.stim_type(count1)) & H.stim_type(count1) == 9		% from file
%         set(H.DUR,'Enable','off');
%         clear dir
%         D = dir([FN.stim_path char(FN.stim_ted(count1))]);
%         XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
%         set(H.DUR,'string',num2str(XStimParams.curr_stimdur));
%         isthereafileflag = 1;    
%     end
% end
% if ~isthereafile_flag
%     set(H.DUR,'Enable','on');
% end



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


% check for visibility of ITD
if XStimParams.ildalone_flag | XStimParams.ABLalone_flag
    set(H.ITD, 'Enable', 'on');
else
    set(H.ITD, 'Enable', 'off');
end

% check for visibility of ILD
if XStimParams.itdalone_flag | XStimParams.ABLalone_flag
    set(H.az, 'Enable', 'on');
    set(H.el, 'Enable', 'on');
else
    set(H.az, 'Enable', 'off');
    set(H.el, 'Enable', 'off');
end

% set ITD, and ABL
XStimParams.curr_ITD = str2num(get(H.ITD, 'string'));
XStimParams.curr_ABL = str2num(get(H.ABL, 'string'));
XStimParams.test_ISI = str2num(get(H.ISI, 'string'));

%check visibility of record datafile
if get(H.recorddata,'Value')
    set(H.recorddata_FN,'enable','on');
else
    set(H.recorddata_FN,'Enable','off');
end

% check out modulation parameters
% XStimParams.mod_depth = str2num(get(H.Twosrc_mod_depth,'String'));
% mod_type = get(H.Twosrc_mod_type,'String');
% mod_num = get(H.Twosrc_mod_type,'Value');
% XStimParams.mod_type = deblank(mod_type(mod_num,:));

% switch XStimParams.mod_type
%     case 'Tone'
% 		set(H.Twosrc_mod_txt,'Visible', 'on');
%         set(H.Twosrc_mod_txt,'String', 'mod freq ');
%         set(H.Twosrc_mod_freq,'Visible','on');
%         set(H.Twosrc_mod_phase,'Visible','on');
%         set(H.Twosrc_mod_txtA,'Visible','on');
%         set(H.Twosrc_mod_depth,'Visible','on');
%         set(H.Twosrc_mod_txtB,'Visible','on');
% 		set(H.Twosrc_mod_txtC,'Visible','off');
%         set(H.Twosrc_mod_pb,'Visible','off');
%         XStimParams.mod_freq = str2num(get(H.Twosrc_mod_freq,'String'));
%         XStimParams.mod_phase = str2num(get(H.Twosrc_mod_phase,'String'));
%     case 'File'
% 		set(H.Twosrc_mod_txt,'Visible', 'off');
%         set(H.Twosrc_mod_pb,'Visible','on');
%         set(H.Twosrc_mod_freq,'Visible','off');
%         set(H.Twosrc_mod_txtA,'Visible','on');
%         set(H.Twosrc_mod_depth,'Visible','on');
%         set(H.Twosrc_mod_phase,'Visible','off');
%         set(H.Twosrc_mod_txtB,'Visible','off');
%         set(H.Twosrc_mod_txtC,'Visible','on');
% 		if get(H.Twosrc_mod_pb,'Value')
%             [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
%             set(H.Twosrc_mod_pb,'Value',0);
%         end
%         set(H.Twosrc_mod_txtC,'String', FN.mod);
%        		
%     case 'LP Noise'
%         set(H.Twosrc_mod_txt,'Visible', 'on');
% 		set(H.Twosrc_mod_txt,'String', 'cutoff freq ');
%         set(H.Twosrc_mod_freq,'Visible','on');
%         set(H.Twosrc_mod_txtA,'Visible','on');
%         set(H.Twosrc_mod_depth,'Visible','on');
%         set(H.Twosrc_mod_phase,'Visible','off');
%         set(H.Twosrc_mod_pb,'Visible','off');
%         set(H.Twosrc_mod_txtB,'Visible','off');
%         set(H.Twosrc_mod_txtC,'Visible','off');
% 		XStimParams.mod_freq = str2num(get(H.Twosrc_mod_freq,'String'));
%         
%     case 'None'
% 		set(H.Twosrc_mod_txt,'Visible', 'on');
%         set(H.Twosrc_mod_txt,'String', 'no mod  ');
%         set(H.Twosrc_mod_freq,'Visible','off');
%         set(H.Twosrc_mod_phase,'Visible','off');
%         set(H.Twosrc_mod_pb,'Visible','off');
%         set(H.Twosrc_mod_txtB,'Visible','off');
%         set(H.Twosrc_mod_txtA,'Visible','off');
%         set(H.Twosrc_mod_depth,'Visible','off');
%         set(H.Twosrc_mod_txtC,'Visible','off');
%     otherwise
% end

% calc Vstrength
if get(H.plotVstr,'value')
    tit = 'Vector Strength Plots';
    prompt = {'Enter datapath', 'Enter filename ','Enter vector of modulation frequencies'};
    DefAns = {FN.data_path,[FN.data '.mat'], '[55 75]'};
    
    Ans = inputdlg(prompt,tit,1,DefAns);
    read_data3_minus_spont(Ans{2},Ans{1}, [], 1,str2num(Ans{3}));
    clear tit prompt DefAns Ans
    set(H.plotVstr,'value',0);
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
update_XStim