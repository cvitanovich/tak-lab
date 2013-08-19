% setInfo_Xstim
% combines callbacks for Xstim
global XStimParams

XStimParams.ephone_flag = get(H.ephoneuseit,'Value');
XStimParams.ap_pos = str2num(get(H.ap_pos,'String'));
XStimParams.ml_pos = str2num(get(H.ml_pos,'String'));
XStimParams.depth = str2num(get(H.depth,'String'));

if ~strcmp(XStimParams.recording_site, get(H.recording_site, 'String'))
    XStimParams.recording_site = get(H.recording_site, 'String');
    tempstr = [num2str(XStimParams.depth) ' um depth'];
    update_diary;
end
XStimParams.testnum = str2num(get(H.testnum, 'String'));

% update XStimParams and FN to birdnum
if XStimParams.bird_number ~= str2num(get(H.controlbirdnumber,'String'))
    XStimParams.bird_number = str2num(get(H.controlbirdnumber,'String'));
    Globals_FN_update;
end

% reset a few params
% leading silence dur
XStimParams.silence_lead = str2num(get(H.Xstim_silence_lead,'string'));
XStimParams.silence_trail = str2num(get(H.Xstim_silence_trail,'string'));
XStimParams.curr_stimdur = str2num(get(H.Xstim_curr_stimdur,'string'));
XStimParams.test_ISI = str2num(get(H.Xstim_test_ISI,'string'));

% get stimulus1 type from XStimParams
stim_type = get(H.stim_type,'String');
Nstims = size(stim_type,1);
for i = 1:Nstims
   if find(strcmp(deblank(stim_type(i,:)), XStimParams.stim_type))
      set(H.stim_type,'Value',i);
   end
end
stim_type = get(H.stim_type,'String');
stim_val = get(H.stim_type,'Value');
if strcmp(deblank(stim_type(stim_val,:)), 'File')
   set(H.stim_filename,'Enable','on','String',FN.stim);
else
   set(H.stim_filename,'Enable','off')
end

% get stimulus2 type from XStimParams
for i = 1:Nstims
   if find(strcmp(deblank(stim_type(i,:)), XStimParams.stim_type2))
      set(H.stim_type2,'Value',i);
   end
end
stim_type = get(H.stim_type2,'String');
stim_val = get(H.stim_type2,'Value');
if strcmp(deblank(stim_type(stim_val,:)), 'File')
   set(H.stim_filename2,'Enable','on','String',FN.stim2)
else
   set(H.stim_filename2,'Enable','off')
end

% get stimulus3 type from XStimParams
for i = 1:Nstims
   if find(strcmp(deblank(stim_type(i,:)), XStimParams.stim_type3))
      set(H.stim_type3,'Value',i);
   end
end
stim_type = get(H.stim_type3,'String');
stim_val = get(H.stim_type3,'Value');
if strcmp(deblank(stim_type(stim_val,:)), 'File')
   set(H.stim_filename3,'Enable','on','String',FN.stim3)
else
   set(H.stim_filename3,'Enable','off')
end

% reset displays in space, altIR and 2-source
if exist1('H.space_stim_type')
    set(H.space_stim_type,'Value',stim_val);
end
if exist1('H.space4_stim_type')
    set(H.space4_stim_type,'Value',stim_val);
end
if exist1('H.Twosrc_stim_type')
    set(H.Twosrc_stim_type,'Value',stim_val);
end

% earphone file
if get(H.ephonepb, 'Value')
    [FN.ephone, FN.ephone_path] = uigetfile([FN.ephone_path '*.imp'],'Select Earphone Filter File');
    if(FN.ephone_path ~= 0)
        set(H.ephonefile,'String',FN.ephone);
        eval(['save ' FN.current_path 'FN_current FN;'])
    end
    FN.HRTFfiletype(5) = testHRTFfiletype(FN.ephone_path, FN.ephone);
end

% earphone file2
if get(H.ephonepb2, 'Value')
    [FN.ephone2, FN.ephone_path] = uigetfile([FN.ephone_path '*.imp'],'Select Earphone Filter File');
    if(FN.ephone_path ~= 0)
        set(H.ephonefile2,'String',FN.ephone2);
        eval(['save ' FN.current_path 'FN_current FN;'])
    end
    FN.HRTFfiletype(6) = testHRTFfiletype(FN.ephone_path, FN.ephone2);
end

if get(H.locpb,'Value')
    [FN.loc, FN.loc_path] = uigetfile([FN.loc_path '*.*'],'Select Location Filter File');
    if(FN.loc_path ~= 0)
        set(H.locfile,'String',[FN.loc_path FN.loc]);
	    eval(['save ' FN.current_path 'XStimParams_current XStimParams']);
        eval(['save ', FN.current_path  'FN_current FN'])
    end
    FN.HRTFfiletype(4) = testHRTFfiletype(FN.loc_path, FN.loc);
end

% increment test number
if get(H.Xstim_inctestnum,'Value')
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    set(H.Xstim_inctestnum,'Value',0);
end
update_dataFN;
set(H.xstim_recorddata_FN, 'String', FN.data);

%SetLocFlag
XStimParams.loc_flag = get(H.locuseit,'Value');
AZ = str2num(get(H.locAZ,'String'));
EL = str2num(get(H.locEL,'String'));
XStimParams.loc_azel = [AZ EL];

% diary entry
if get(H.diaryEntry,'value')
    set(H.diaryEntry,'value',0);
    tempstr = inputdlg('Enter text for diary','Diary Entry');
    tempstr = tempstr{1};
    update_diary;
end

% S2close
if get(H.S2closepb,'value')
    global_exit
end

% ClearStim Directory
if get(H.ClearStimspb,'value')
    eval(['delete ' FN.temp_stim_path '*.*;']);    
end

% check out test type
test_type = get(H.test_type,'String');
stim_val = get(H.test_type,'Value');
XStimParams.test_type = deblank(test_type(stim_val,:));


if(exist1('H.ablfig') | ...
exist1('H.searchfig') | ...
exist1('H.spikefig') | ...
exist1('H.itdfig') | ...
exist1('H.itd2fig') | ...
exist1('H.itd_decorrfig') | ...
exist1('H.freqtestfig') | ...
exist1('H.ildfreqfig') | ...
exist1('H.spacefig') | ...
exist1('H.spacefig3') | ...
exist1('H.space4fig') | ...
exist1('H.Two_sourcefig') | ...
exist1('H.altIRfig') | ...
exist1('H.Delayfig') | ...
exist1('H.multi_sourcefig') | ...
exist1('H.AM_fig') | ...
exist1('H.PE_envfig'))

else


switch stim_val
    case 1
    case 2     % Search for cells
        searchcells;
    case 3     % ABL Test
        abltest;
    case 4     % ITD Test
        itdtest;
    case 5     % ILD Test
        ildtest;
    case 6     % FREQ Test
        freqtest;
    case 7     % Dean & McAlpine like test co-localized
        mcSpace;
    case 8      % Dean & McAlpine like test 2 sources
        mc2Source;
    case 9     % Space no dbl buff, *.std HRIRs, all stims saved to disk (already convolved with HRIR) before playout
        space;
    case 10     % space2 no dbl buff, *.std HRIRs, does NOT write files to disk
        space2;
    case 11 	 % space3 with dbl buff, uses *.eq HRIRs, L&R stim written to disk before playout
        space3;
    case 12		 % basic 2 source test
        Two_Source;
    case 13		 % alterred IR (Hanna's)
        alt_IR_XXX;  
    case 14
        AM;
    case 15         % allows filename for each ear
        itdtest2;
    case 16         % tests itds with different decorrelations (from files)
        itd_decorr;
    case 17
        composite;  % plays out saved stims in temp_stim_dir
    case 18
        BMLD;
    case 19
        LRsounds;       % plays files
    case 20
        Adapt1;         % adaptation
    case 21
        Adapt2;         % adaptation for two sources
    case 22
        space4;         % plays Nagel & Doupe-like stimuli
    case 23
        McSpace2;       % plays two different Dean/McAlpine regimes
    case 24
        Space_FileperRep;   % plays nReps of each file chosen
    case 25
        LongSAMs;           % plays SAM noise via play2_record2b_SAMB.dll
    case 26
        MaskedSpace;        % uses DLL to play masker then probe
    case 27
        PE_env;             % precedence effect with envelopes from gammatoned noise
    case 28
        Space_RIR;          % room impulse responses
    case 29
        AdaptedSpace;       % adaptor-probe presentation with fully-cued probe locations
    otherwise
        disp('Not implemented')
        set(H.test_type,'Value',1);
end
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
%eval(['save ' FN.current_path 'FN_current FN;'])
