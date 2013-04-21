% SetInfo_altIR
% combines all the other 'set...' scripts into one

%check ILDAloneFlag
% XStimParams.space_flag == 1 for ILDalone
%								 == 0 for fully cued

% uses .STD HRIRs, therefore make sure earphone filters are also selected
   XStimParams.ephone_flag = 1;
   set(H.ephoneuseit,'Value',1);

      
if get(H.altIR_buildOnly,'value') | exist1('H.composite_fig')
    XStimParams.buildOnly = 1;
else
    XStimParams.buildOnly = 0;
end   
set(H.altIR_buildOnly,'value',XStimParams.buildOnly);

if XStimParams.space_flag == get(H.altIR_ildaloneflag,'value') | ...
   XStimParams.space_flag ~= get(H.altIR_spaceflag,'value')   % value has changed
   if XStimParams.space_flag == 0		% was ILDalone, now fully cued
		set(H.altIR_spaceflag,'Value',1);
        set(H.altIR_ildaloneflag,'Value',0);
        XStimParams.ildalone_flag = 0;
   else									% was fully cued, now ILDalone
		set(H.altIR_spaceflag,'Value',0);
        set(H.altIR_ildaloneflag,'Value',1);
        XStimParams.ildalone_flag = 1;
   end
   XStimParams.space_flag = ~XStimParams.space_flag;
end

% control visibility of file picking for ILA or fully-cued
if get(H.altIR_spaceflag,'Value')
   set(H.altIR_spacefilehdr,'Enable','on');
   set(H.altIR_spacefile,'Enable','on');
   set(H.altIR_spacefilepb,'Enable','on');   
   set(H.altIR_ildalonefilehdr,'Enable','off');
   set(H.altIR_ildalonefile,'Enable','off');
   set(H.altIR_ildalonepb,'Enable','off');  
   set(H.altIR_ildaloneflag,'Value',0);
   set(H.altIR_ITDchk,'Value',0);
else
   set(H.altIR_spacefilehdr,'Enable','off');
   set(H.altIR_spacefile,'Enable','off');
   set(H.altIR_spacefilepb,'Enable','off');
   set(H.altIR_ildalonefilehdr,'Enable','on');
   set(H.altIR_ildalonefile,'Enable','on');
   set(H.altIR_ildalonepb,'Enable','on');
   set(H.altIR_ildaloneflag,'value',1);
   set(H.altIR_ITDchk,'Value',1);
end

% check for IA pushbutton (file-picker)
if get(H.altIR_ildalonepb,'Value')
   [FN.ildalone_std,FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select ILD-alone HRTF File');
	if(FN.ILA_path ~= 0)
   	set(H.altIR_ildalonefile,'String',FN.ildalone_std);
   end
   set(H.altIR_ildalonepb,'Value',0);
    FN.HRTFfiletype(2) = testHRTFfiletype(FN.ILA_path, FN.ildalone_std);
end
   
% check for fully-cued pushbutton (file-picker)
if get(H.altIR_spacefilepb,'Value')
    [FN.space_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File (*.std preferred)');
  	if(FN.space_path ~= 0)
   	set(H.altIR_spacefile,'String',FN.space_std);
  	end
   set(H.altIR_spacefilepb,'Value',0);
    FN.HRTFfiletype(1) = testHRTFfiletype(FN.space_path, FN.space_std);
end

% check for saving or using saved altered HRTFs or ALTerations
if get(H.altIR_saveHRTFflag,'Value')  
    XStimParams.use_saveHRTFflag = 0;
    set(H.altIR_use_saveHRTFflag,'Value',0);
    XStimParams.use_saveALTflag = 0;
    set(H.altIR_use_saveALTflag,'Value',0);
end
if get(H.altIR_use_saveHRTFflag,'Value')
    XStimParams.use_saveHRTFflag = 1;
    set(H.altIR_saveHRTFflag,'Value',0);
    XStimParams.use_saveALTflag = 0;
    set(H.altIR_use_saveALTflag,'Value',0);
end
if get(H.altIR_use_saveALTflag,'Value')
    XStimParams.use_saveALTflag = 1;
    set(H.altIR_saveHRTFflag,'Value',0);
    XStimParams.use_saveHRTFflag = 0;
    set(H.altIR_use_saveHRTFflag,'Value',0);
end

if get(H.altIR_saveHRTFflag,'Value')
    while isempty(FN.saveHRTF) | exist([FN.HRTF_path FN.saveHRTF{1} '.mat'],'file')
       FN.saveHRTF = inputdlg('Enter New Filename to SAVE: ');
       set(H.altIR_saveHRTF_FN,'String', FN.saveHRTF);
    end
end
if get(H.altIR_use_saveHRTFflag,'Value') | get(H.altIR_use_saveALTflag,'Value')
    while isempty(FN.saveHRTF) | ~exist([FN.HRTF_path FN.saveHRTF{1} '.mat'],'file')
       FN.saveHRTF = inputdlg('Enter New Filename to USE: ');
       set(H.altIR_saveHRTF_FN,'String', FN.saveHRTF);
    end
end

if ~strcmp(get(H.altIR_saveHRTF_FN,'String'),FN.saveHRTF)
    FN.saveHRTF = get(H.altIR_saveHRTF_FN,'String');
end
set(H.altIR_saveHRTF_FN,'String', FN.saveHRTF);

% set numreps
XStimParams.numreps = str2num(get(H.altIR_numreps,'String'));

%Set space Picker
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
if(exist1('H.altIRfig') & get(H.altIR_inc_testnum,'Value'))
	XStimParams.testnum = XStimParams.testnum +1;
	set(H.testnum, 'String', num2str(XStimParams.testnum))
   update_dataFN;
   set(H.altIR_inc_testnum,'Value',0);
end

% update display of stim filename
if get(H.altIR_stimFNpb, 'Value')
    Stimulus_type1;
    set(H.altIR_stimFNpb, 'Value',0);
    set(H.stim_type,'Value',9);
end
if get(H.stim_type,'Value') == 9								% from file
   set(H.altIR_stimFN,'String',FN.stim);
   set(H.stim_filename,'String',FN.stim);
   set(H.stim_filename,'Enable','on');
   set(H.altIR_stimFN,'Enable','on');
else
   set(H.stim_filename,'Enable','off');
   set(H.altIR_stimFN,'Enable','off');
end

% update display of recorddata FN
set(H.altIR_recorddata_FN,'String', FN.data);

% LEFT & RIGHT
XStimParams.alter_left = get(H.alter_left,'Value');
XStimParams.alter_right = get(H.alter_right,'Value');

%MAGNITUDE
XStimParams.alterMag = get(H.alterMag,'Value');
XStimParams.alterMag_dB = str2num(get(H.alterMag_dB,'String'));
%PHASE
XStimParams.alterPhase = get(H.alterPhase,'Value');
XStimParams.alterPhase_rad = str2num(get(H.alterPhase_rad,'String'));

% get octlists
XStimParams.Twelths_to_combine = str2num(get(H.altIR_Twelths_to_combine,'string'));
if ~ismember(XStimParams.Twelths_to_combine,[1 2 3 4 6])
   uiwait(helpdlg('Twelvths_to_combine must be one of [1 2 3 4 6]','PLEASE CHECK #Twelvths'));
end
for i = 1:n_cF-1
    XStimParams.Mag_octlistBIN(i) = get(H.altIR_Mag_octlist(i),'value');
    XStimParams.Phase_octlistBIN(i) = get(H.altIR_Phase_octlist(i),'value');
end

if 0
XStimParams.alter_setbands = get(H.alter_setbands,'value');
if ~isempty(get(H.alterMag_Minbandlimits(1),'String'))   
   XStimParams.alterMag_Minbandlimits(1) = str2num(get(H.alterMag_Minbandlimits(1),'String'));
end
    
if ~isempty(get(H.alterMag_Maxbandlimits(1),'String'))   
   XStimParams.alterMag_Maxbandlimits(1) = str2num(get(H.alterMag_Maxbandlimits(1),'String'));
end
    
if ~isempty(get(H.alterPhase_Minbandlimits(1),'String'))   
   XStimParams.alterPhase_Minbandlimits(1) = str2num(get(H.alterPhase_Minbandlimits(1),'String'));
end
    
if ~isempty(get(H.alterPhase_Maxbandlimits(1),'String'))   
   XStimParams.alterPhase_Maxbandlimits(1) = str2num(get(H.alterPhase_Maxbandlimits(1),'String'));
end
end

if ~isempty(XStimParams.locations)
    locations = XStimParams.locations;
    eval(['save ' FN.current_path 'Locations_current locations;'])
end
%clear locations
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim