%Stimulus_type3: Callback for Stimulus Type2 popupmenu

stim_type = get(H.stim_type3,'String');
stim_val = get(H.stim_type3,'Value');
XStimParams.stim_type3 = deblank(stim_type(stim_val,:));
clear stim_type

if(strcmp('File',deblank(XStimParams.stim_type2)))
   [stim_from_file3, FN.stim3, FN.stim_path3] = Get_Stimulus_File(XStimParams, FN.stim_path3);
   %XStimParams.stimfile3 = [FN.stim_path3 FN.stim3];
end

% show filename on Xstim control panel
if get(H.stim_type3,'Value')==9
   set(H.stim_filename3,'Enable','on',...
      'String', FN.stim3);
else
  	set(H.stim_filename3,'Enable','off');   
end
eval(['save ' FN.current_path 'XStimParams_current XStimParams']);
