%Stimulus_type1: Callback for Stimulus Type1 popupmenu

stim_type = get(H.stim_type,'String');
stim_val = get(H.stim_type,'Value');
XStimParams.stim_type = deblank(stim_type(stim_val,:));
clear stim_type

if(strcmp('File',deblank(XStimParams.stim_type)))
   [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
   %XStimParams.stimfile = [FN.stim_path FN.stim];
end

% reset displays in space and 2-source
if exist1('H.space_stim_type')
    set(H.space_stim_type,'Value',stim_val);
end
if exist1('H.Twosrc_stim_type')
    set(H.Twosrc_stim_type,'Value',stim_val);
end


% show filename on Xstim control panel
if get(H.stim_type,'Value')==9
   set(H.stim_filename,'Enable','on',...
      'String', FN.stim);
else
  	set(H.stim_filename,'Enable','off');   
end
eval(['save ' FN.current_path 'XStimParams_current XStimParams']);

if (exist1('H.spacefig3'))
   setInfo_space3;
end
