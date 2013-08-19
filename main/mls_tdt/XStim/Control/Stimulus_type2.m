%Stimulus_type2: Callback for Stimulus Type2 popupmenu

stim_type = get(H.stim_type2,'String');
stim_val = get(H.stim_type2,'Value');
XStimParams.stim_type2 = deblank(stim_type(stim_val,:));
clear stim_type

if(strcmp('File',deblank(XStimParams.stim_type2)))
   [stim_from_file2, FN.stim2, FN.stim_path2] = Get_Stimulus_File(XStimParams, FN.stim_path2);
   %XStimParams.stimfile2 = [FN.stim_path2 FN.stim2];
end

% show filename on Xstim control panel
if get(H.stim_type2,'Value')==9
   set(H.stim_filename2,'Enable','on',...
      'String', FN.stim2);
else
  	set(H.stim_filename2,'Enable','off');   
end
eval(['save ' FN.current_path 'XStimParams_current XStimParams']);
