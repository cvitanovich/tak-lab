%InitXStimParams: Callback to pushbutton H.initsessgo

XStimParams.testnum = str2num(get(H.testnum_init,'String'));
XStimParams.recording_site = get(H.recording_site_init,'string');
XStimParams.datadir = FN.data_path;

% reset bird_number, ephone use and data filename
if(get(H.paramdefaults,'Value'))			% use defaults
	eval(['load ', FN.current_path  'XStimParams_default'])
   XStimParams.bird_number = str2num(get(H.birdnuminit,'String'));
   XStimParams.ephonefile = get(H.ephonefile,'String');
   update_dataFN;
	eval(['save ' FN.current_path 'XStimParams_default XStimParams']);
else								% use current values
   XStimParams.bird_number = str2num(get(H.birdnuminit,'String'));
   update_dataFN;
end

XStimParams.datafile = FN.data;

eval(['save ' FN.current_path 'XStimParams_current XStimParams']);

set(H.controlbirdnumber,'String',num2str(XStimParams.bird_number));
set(H.ap_pos,'String',num2str(XStimParams.ap_pos));
set(H.ml_pos,'String',num2str(XStimParams.ml_pos));
set(H.depth,'String',num2str(XStimParams.depth));
set(H.ephoneuseit,'Value',XStimParams.ephone_flag);
set(H.locuseit,'Value',XStimParams.loc_flag);
set(H.testnum,'String',num2str(XStimParams.testnum));
set(H.recording_site, 'String',XStimParams.recording_site);

testtype = get(H.test_type,'String');
for s=1:size(testtype,1)
   if(strcmp(deblank(testtype(s,:)),XStimParams.test_type))
      set(H.test_type,'Value',s);
   end
end

stim_type = get(H.stim_type,'String');
for s=1:size(stim_type,1)
   if(strcmp(deblank(stim_type(s,:)),XStimParams.stim_type))
      set(H.stim_type,'Value',s);
   end
end

close(H.initsessfig)
H.initsessfig = [];
