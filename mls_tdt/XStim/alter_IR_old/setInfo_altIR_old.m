% SetInfo_altIR
% combines all the other 'set...' scripts into one

%check ILDAloneFlag
% XStimParams.space_flag == 1 for ILDalone
%								 == 0 for fully cued

if XStimParams.space_flag == get(H.ildaloneflag,'value') | ...
   XStimParams.space_flag ~= get(H.spaceflag,'value')   % value has changed
   if XStimParams.space_flag == 0		% was ILDalone, now fully cued
		set(H.spaceflag,'Value',1);
   else									% was fully cued, now ILDalone
		set(H.spaceflag,'Value',0);
   end
   XStimParams.space_flag = ~XStimParams.space_flag;
end

% control visibility of file picking for ILA or fully-cued
if get(H.spaceflag,'Value')
   set(H.spacefilehdr,'Enable','on');
   set(H.spacefile,'Enable','on');
   set(H.spacefilepb,'Enable','on');   
   set(H.ildalonefilehdr,'Enable','off');
   set(H.ildalonefile,'Enable','off');
   set(H.ildalonepb,'Enable','off');  
   set(H.ildaloneflag,'Value',0);
   set(H.ITDchk,'Value',0);
%Make sure earphone filters are not selected
   XStimParams.ephone_flag = 0;
   set(H.ephoneuseit,'Value',0);
else
   set(H.spacefilehdr,'Enable','off');
   set(H.spacefile,'Enable','off');
   set(H.spacefilepb,'Enable','off');
   set(H.ildalonefilehdr,'Enable','on');
   set(H.ildalonefile,'Enable','on');
   set(H.ildalonepb,'Enable','on');
   set(H.ildaloneflag,'value',1);
   set(H.ITDchk,'Value',1);
%Make sure earphone filters are also selected
   XStimParams.ephone_flag = 1;
   set(H.ephoneuseit,'Value',1);
end

% check for IA pushbutton (file-picker)
if get(H.ildalonepb,'Value')
   [FN.ila,FN.ILA_path] = uigetfile('*.*','Select ILD-alone HRTF File');
	if(FN.ILA_path ~= 0)
   	set(H.ildalonefile,'String',FN.ila);
   end
   set(H.ildalonepb,'Value',0);
end
   
% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
  	[FN.space,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File');
  	if(FN.space_path ~= 0)
   	set(H.spacefile,'String',FN.space);
  	end
   set(H.spacefilepb,'Value',0);
end

% set numreps
XStimParams.numreps = get(H.numreps,'String');

%Set space PickerFlag
if get(H.pickerflag,'Value')
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
   set(H.locscriptfilehdr,'Enable','off');
   set(H.locscriptfile,'Enable','off');
   set(H.locscriptfilepb,'Enable','off');
else
   XStimParams.picker_flag = 0;
   if exist1('H.pickerfig') 
      close(H.pickerfig); 
      H.pickerfig = [];
      XStimParams.locations = [];
   end
   set(H.locscriptfilehdr,'Enable','on');
   set(H.locscriptfile,'Enable','on');
   set(H.locscriptfilepb,'Enable','on');
end

%draw_loc_scrpt

% increment test number
if(exist1('H.altIRfig') & get(H.inc_testnum,'Value'))
	XStimParams.testnum = XStimParams.testnum +1;
	set(H.testnum, 'String', num2str(XStimParams.testnum))
   update_dataFN;
   set(H.inc_testnum,'Value',0);
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

for i = 1:3     % fill XStimParams. ... by reading strings from controlfig
    if ~isempty(get(H.alterMag_Minbandlimits(i),'String'))   
        XStimParams.alterMag_Minbandlimits(i) = str2num(get(H.alterMag_Minbandlimits(i),'String'));
    elseif length(XStimParams.alterMag_Minbandlimits) >= i
        XStimParams.alterMag_Minbandlimits(i) = [];
    end
    
    if ~isempty(get(H.alterMag_Maxbandlimits(i),'String'))   
        XStimParams.alterMag_Maxbandlimits(i) = str2num(get(H.alterMag_Maxbandlimits(i),'String'));
    elseif length(XStimParams.alterMag_Maxbandlimits) >= i
        XStimParams.alterMag_Maxbandlimits(i) = [];
    end
    
    if ~isempty(get(H.alterPhase_Minbandlimits(i),'String'))   
        XStimParams.alterPhase_Minbandlimits(i) = str2num(get(H.alterPhase_Minbandlimits(i),'String'));
    elseif length(XStimParams.alterPhase_Minbandlimits) >= i
        XStimParams.alterPhase_Minbandlimits(i) = [];
    end
    
    if ~isempty(get(H.alterPhase_Maxbandlimits(i),'String'))   
        XStimParams.alterPhase_Maxbandlimits(i) = str2num(get(H.alterPhase_Maxbandlimits(i),'String'));
    elseif length(XStimParams.alterPhase_Maxbandlimits) >= i
        XStimParams.alterPhase_Maxbandlimits(i) = [];
    end
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])

