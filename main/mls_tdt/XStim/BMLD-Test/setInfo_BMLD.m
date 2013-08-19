% SetInfo_BMLD
% combines all the other 'set...' scripts into one

%check ILDAloneFlag
% XStimParams.space_flag == 1 for ILDalone
%						 == 0 for fully cued

XStimParams.ildalone_flag = 0;
XStimParams.itdalone_flag = 0;
XStimParams.ABLalone_flag = 0;
if get(H.spaceflag,'Value')
    XStimParams.space_flag = 1;
else
    XStimParams.space_flag = 0;
end
if get(H.each_loc,'Value')
    XStimParams.each_loc = 1;
else
    XStimParams.each_loc = 0;
end

% ITD/locations visible or not
if get(H.spaceflag,'Value')
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');
    for iSrc = 1:2
        set(H.ITD_txt(iSrc),'Enable','off');
        set(H.ITD(iSrc),'Enable','off');
        set(H.loc1_txt(iSrc),'Enable','on');
        set(H.loc1(iSrc),'Enable','on');
        set(H.loc2_txt(iSrc),'Enable','on');
        set(H.loc2(iSrc),'Enable','on');
    end
    XStimParams.locations(1) = str2num(get(H.loc1(1),'string'));
    XStimParams.locations(2) = str2num(get(H.loc1(2),'string'));
    XStimParams.locations1(1) = str2num(get(H.loc2(1),'string'));
    XStimParams.locations1(2) = str2num(get(H.loc2(2),'string'));
    
else
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');   
    for iSrc = 1:2
        set(H.ITD_txt(iSrc),'Enable','on');
        set(H.ITD(iSrc),'Enable','on');
        set(H.loc1_txt(iSrc),'Enable','off');
        set(H.loc1(iSrc),'Enable','off');
        set(H.loc2_txt(iSrc),'Enable','off');
        set(H.loc2(iSrc),'Enable','off');
    end
    XStimParams.curr_ITD = str2num(get(H.ITD(1),'string'));
end
    

% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
    [FN.space_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space_std);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
end

% set numreps
XStimParams.numreps = str2num(get(H.numreps,'String'));

% increment test number
if(exist1('H.BMLDfig') & get(H.inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.inc_testnum,'Value',0);
end

% update stim_type1
stim_type = get(H.BMLD_stim_type,'String');
stim_val = get(H.BMLD_stim_type,'Value');

if ~strcmp(XStimParams.stim_type, deblank(stim_type(stim_val,:))) | ...
        get(H.BMLD_stim_filenamepb,'Value')		% stim_type has changed or pb pushed
    XStimParams.stim_type = deblank(stim_type(stim_val,:));
    set(H.stim_type,'Value',stim_val);
    if(strcmp('File',deblank(XStimParams.stim_type)))
        [stim_from_file, FN.stim, FN.stim_path] = Get_Stimulus_File(XStimParams, FN.stim_path);
        set(H.BMLD_stim_filename,'String', FN.stim);
        set(H.stim_filename,'String', FN.stim);
    end
    %setinfo_Xstim;
    set(H.BMLD_stim_filenamepb,'Value',0);
end

%control visibility of FN vs freq
switch deblank(XStimParams.stim_type)
    case 'File'
        set(H.BMLD_stim_filenamepb,'Visible','on');
        set(H.BMLD_stim_filename,'Visible','on');
        set(H.BMLD_freq_txt,'Visible','off');
        set(H.BMLD_freq,'Visible','off');    
    case {'Gammatone', 'Tone'}
        set(H.BMLD_freq_txt,'Visible','on');
        set(H.BMLD_freq,'Visible','on');    
        set(H.BMLD_stim_filenamepb,'Visible','off');
        set(H.BMLD_stim_filename,'Visible','off');
    otherwise
        set(H.BMLD_stim_filenamepb,'Visible','off');
        set(H.BMLD_stim_filename,'Visible','off');
        set(H.BMLD_freq_txt,'Visible','off');
        set(H.BMLD_freq,'Visible','off');    
end

% update stim_type2
stim_type = get(H.BMLD_stim_type2,'String');
stim_val = get(H.BMLD_stim_type2,'Value');

if ~strcmp(XStimParams.stim_type2,deblank(stim_type(stim_val,:))) | ...
        get(H.BMLD_stim_filename2pb,'Value')		% stim_type2 has changed
    XStimParams.stim_type2 = deblank(stim_type(stim_val,:));
    set(H.stim_type2,'Value',stim_val);
    if(strcmp('File',deblank(XStimParams.stim_type2)))
        [stim_from_file2, FN.stim2, FN.stim_path2] = Get_Stimulus_File(XStimParams, FN.stim_path2);
        set(H.BMLD_stim_filename2,'String', FN.stim2);
        set(H.stim_filename2,'String', FN.stim2);
    end
    set(H.BMLD_stim_filename2pb,'Value',0);
end
%control visibility of FN vs freq
switch deblank(XStimParams.stim_type2)
    case 'File'
        set(H.BMLD_stim_filename2,'Enable','on');
        set(H.stim_filename2,'Enable','on');
        set(H.BMLD_stim_filename2pb,'Visible','on');
        set(H.BMLD_stim_filename2,'Visible','on');
        set(H.BMLD_freq2_txt,'Visible','off');
        set(H.BMLD_freq2,'Visible','off');    
    case {'Gammatone', 'Tone'}
        set(H.BMLD_freq2_txt,'Visible','on');
        set(H.BMLD_freq2,'Visible','on');    
        set(H.BMLD_stim_filename2pb,'Visible','off');
        set(H.BMLD_stim_filename2,'Visible','off');
    otherwise
        set(H.BMLD_stim_filename2pb,'Visible','off');
        set(H.BMLD_stim_filename2,'Visible','off');
        set(H.BMLD_freq2_txt,'Visible','off');
        set(H.BMLD_freq2,'Visible','off');    
end

% set dur based on stim1
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
if get(H.stim_type,'Value') == 9		% from file
    set(H.DUR,'Enable','off');
    clear dir
    D = dir([FN.stim_path FN.stim]);
    XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
    set(H.DUR,'string',num2str(XStimParams.curr_stimdur));
else
    set(H.DUR,'Enable','on');
end

% bandwidths and bandlimits
% 1st stim
if ~strcmp(XStimParams.stim_type, 'File') 	% if not from file
    % if changed, re-set bandwidths
    if XStimParams.bandwidth(1) ~= str2num(get(H.bandwidth1,'String'))
        XStimParams.bandwidth(1) = str2num(get(H.bandwidth1,'String'));
    end
    % if changed, re-set bandlimits
    if (XStimParams.bandlimits(1) ~= str2num(get(H.bandlimits1,'String')) | ...
            XStimParams.bandlimits(2) ~= str2num(get(H.bandlimits2,'String')));  
        XStimParams.bandlimits(1) = str2num(get(H.bandlimits1,'String'));
        XStimParams.bandlimits(2) = str2num(get(H.bandlimits2,'String'));
    end
    set(H.bandwidth1,'Enable','on');
    set(H.bandlimits1,'Enable','on');
    set(H.bandlimits2,'Enable','on');
else
    set(H.bandwidth1,'Enable','off');
    set(H.bandlimits1,'Enable','off');
    set(H.bandlimits2,'Enable','off');
end

% 2nd stim
if ~strcmp(XStimParams.stim_type2, 'File') 	% if not from file
    % if changed, re-set bandwidths
    if XStimParams.bandwidth(2) ~= str2num(get(H.bandwidth2,'String'))
        XStimParams.bandwidth(2) = str2num(get(H.bandwidth2,'String'));
    end
    % if changed, re-set bandlimits
    if (XStimParams.bandlimits(3) ~= str2num(get(H.bandlimits3,'String')) | ...
            XStimParams.bandlimits(4) ~= str2num(get(H.bandlimits4,'String')));  
        XStimParams.bandlimits(3) = str2num(get(H.bandlimits3,'String'));
        XStimParams.bandlimits(4) = str2num(get(H.bandlimits4,'String'));
    end
    set(H.bandwidth2,'Enable','on');
    set(H.bandlimits3,'Enable','on');
    set(H.bandlimits4,'Enable','on');
else
    set(H.bandwidth2,'Enable','off');
    set(H.bandlimits3,'Enable','off');
    set(H.bandlimits4,'Enable','off');
end

% set ITD, and ABL
XStimParams.curr_ABL = str2num(get(H.ABL,'string'));
XStimParams.test_ISI = str2num(get(H.ISI,'string'));

% set RELATIVE ABLS for 2nd source
if ~get(H.aritrary_abls,'value')
    low_abl  = str2num(get(H.lowabl,'String'));             % should be a negative number
    high_abl = str2num(get(H.highabl,'String'));            % should be a positive number
    if(high_abl + XStimParams.curr_ABL > -20)
        high_abl = -XStimParams.curr_ABL-20;
        set(H.highabl,'String',num2str(high_abl));
    end
    XStimParams.numabls = str2num(get(H.numabls,'String'));
    XStimParams.hiabl = high_abl;
    XStimParams.loabl = low_abl;
    XStimParams.numabls = str2num(get(H.numabls,'String'));   
    set(H.stepabl,'String',num2str(round((high_abl-low_abl)/(XStimParams.numabls-1))));
    XStimParams.abls = round(XStimParams.loabl:...
        (XStimParams.hiabl - XStimParams.loabl)/(XStimParams.numabls-1):...
        XStimParams.hiabl);
    set(H.ablvector,'string',num2str(XStimParams.abls));
else
    XStimParams.abls = str2num(get(H.ablvector,'string'));    
    XStimParams.numabls = length(XStimParams.abls);
    set(H.lowabl,'String','');
    set(H.highabl,'String','');
    set(H.stepabl,'String','');
end

set(H.numabls,'string',num2str(XStimParams.numabls));

%check visibility of record datafile
if get(H.recorddata,'Value')
    set(H.recorddata_FN,'enable','on');
    set(H.recorddata_FN,'string',FN.data);
else
    set(H.recorddata_FN,'Enable','off');
end

% check out modulation parameters for stim1
XStimParams.mod_depth(1) = str2num(get(H.BMLD_mod_depth,'String'));
mod_type = get(H.BMLD_mod_type,'String');
mod_num = get(H.BMLD_mod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Tone'
        set(H.BMLD_mod_txt,'String', 'Mod Freq:');
        set(H.BMLD_mod_freq,'Visible','on');
        set(H.BMLD_mod_phase,'Visible','on');
        set(H.BMLD_mod_txtA,'Visible','on');
        set(H.BMLD_mod_depth,'Visible','on');
        set(H.BMLD_mod_txtB,'Visible','on');
        set(H.BMLD_mod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.BMLD_mod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.BMLD_mod_phase,'String'));
    case 'File'
        set(H.BMLD_mod_pb,'Visible','on');
        set(H.BMLD_mod_freq,'Visible','off');
        set(H.BMLD_mod_txtA,'Visible','on');
        set(H.BMLD_mod_depth,'Visible','on');
        set(H.BMLD_mod_phase,'Visible','off');
        set(H.BMLD_mod_txtB,'Visible','off');
        if get(H.BMLD_mod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.BMLD_mod_pb,'Value',0);
        end
        set(H.BMLD_mod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.BMLD_mod_txt,'String', 'CutOff Freq:');
        set(H.BMLD_mod_freq,'Visible','on');
        set(H.BMLD_mod_txtA,'Visible','on');
        set(H.BMLD_mod_depth,'Visible','on');
        set(H.BMLD_mod_phase,'Visible','off');
        set(H.BMLD_mod_pb,'Visible','off');
        set(H.BMLD_mod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.BMLD_mod_freq,'String'));
        
    case 'None'
        set(H.BMLD_mod_txt,'String', 'no mod  ');
        set(H.BMLD_mod_freq,'Visible','off');
        set(H.BMLD_mod_phase,'Visible','off');
        set(H.BMLD_mod_pb,'Visible','off');
        set(H.BMLD_mod_txtB,'Visible','off');
        set(H.BMLD_mod_txtA,'Visible','off');
        set(H.BMLD_mod_depth,'Visible','off');
    otherwise
end

% check out modulation parameters for stim2
XStimParams.mod_depth(2) = str2num(get(H.BMLD_mod_depth2,'String'));
mod_type = get(H.BMLD_mod_type2,'String');
mod_num = get(H.BMLD_mod_type2,'Value');
XStimParams.mod_type2 = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type2
    case 'Tone'
        set(H.BMLD_mod_txt2,'String', 'Mod Freq:');
        set(H.BMLD_mod_freq2,'Visible','on');
        set(H.BMLD_mod_pb2,'Visible','off');
        set(H.BMLD_mod_txtA2,'Visible','on');
        set(H.BMLD_mod_depth2,'Visible','on');
        set(H.BMLD_mod_phase2,'Visible','on');
        set(H.BMLD_mod_txtB2,'Visible','on');
        XStimParams.mod_freq(2) = str2num(get(H.BMLD_mod_freq2,'String'));
        XStimParams.mod_phase(2) = str2num(get(H.BMLD_mod_phase2,'String'));
    case 'File'
        set(H.BMLD_mod_pb2,'Visible','on');
        set(H.BMLD_mod_freq2,'Visible','off');
        set(H.BMLD_mod_txtA2,'Visible','on');
        set(H.BMLD_mod_depth2,'Visible','on');
        set(H.BMLD_mod_phase2,'Visible','off');
        set(H.BMLD_mod_txtB2,'Visible','off');
        if get(H.BMLD_mod_pb2,'Value')
            [mod_from_file2, FN.mod2, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.BMLD_mod_pb2,'Value',0);
        end
        set(H.BMLD_mod_txt2,'String', FN.mod2);
        
    case 'LP Noise'
        set(H.BMLD_mod_txt2,'String', 'CutOff Freq:');
        set(H.BMLD_mod_freq2,'Visible','on');
        set(H.BMLD_mod_txtA2,'Visible','on');
        set(H.BMLD_mod_depth2,'Visible','on');
        set(H.BMLD_mod_pb2,'Visible','off');
        set(H.BMLD_mod_phase2,'Visible','off');
        set(H.BMLD_mod_txtB2,'Visible','off');
        XStimParams.mod_freq(2) = str2num(get(H.BMLD_mod_freq2,'String'));
        
    case 'None'
        set(H.BMLD_mod_txt2,'String', 'no mod  ');
        set(H.BMLD_mod_freq2,'Visible','off');
        set(H.BMLD_mod_phase2,'Visible','off');
        set(H.BMLD_mod_pb2,'Visible','off');
        set(H.BMLD_mod_txtB2,'Visible','off');
        set(H.BMLD_mod_txtA2,'Visible','off');
        set(H.BMLD_mod_depth2,'Visible','off');   
    otherwise
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])
update_XStim