function [] = Engage_Space3()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_Space3
%*******************************************************************************
%	The Space Test operation Using DSPs to convolve space
%*******************************************************************************
% uses *.eq hrtf files, does not use ephonefilts

disp('This is a SPACE3 TEST and uses DSPs to convolve HRIRs.');

% set params
TDT.itd_scale = .8;
TDT.hrtf_scale = .99;     %.4;
TDT.buf_pts = 300;         %512*8;
DUR_silence = XStimParams.silence_lead; 				%ms

% eliminate files in stim directory
stimuli_dir = FN.temp_stim_path;
fclose all;
if ~get(H.buildOnly,'value')
    eval(['delete ' stimuli_dir '*.*;']);
end

% reset non-allowed stim_type to BroadBand
if isempty(find([1 2 6 8 9] == get(H.space3_stim_type,'Value')))
    set(H.space3_stim_type,'Value',8);
    disp('Unsupported stimulus type for space test; reset to BROADBAND');
end

% get locations from spacepicker
XStimParams.locations = GUI.locations1';
if ~length(XStimParams.locations)
    set(H.pickerfig,'Color', [.1 .2 .8]);
    set(H.picker_error,'visible','on');
    pause;
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end

%Put parameters into XStimParams
test_val = get(H.test_type,'Value');
test_type = get(H.test_type,'String');
XStimParams.test_type = deblank(test_type(test_val,:));
clear test_type
XStimParams.curr_freq = round(str2num(get(H.space3_freq,'String')));
XStimParams.curr_ITD = round(str2num(get(H.space3_ITD,'String')));
XStimParams.curr_ABL = str2num(get(H.space3_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.space3_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space3_ISI,'String'));
XStimParams.numreps = str2num(get(H.space3_numreps,'String'));
XStimParams.az = str2num(get(H.az,'String'));
XStimParams.el = str2num(get(H.el,'String'));
XStimParams.reset_flag = 0;

if(abs(XStimParams.curr_ITD) > 250) disp('ITD is out of range');    return; end

if XStimParams.space_flag
    str1 = 'FC';
    if isempty(FN.space_eq) | ~isempty(strfind(FN.space_eq,'std'))      % FN not yet picked
        [FN.space_eq,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.eq preferred)');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space_eq]);
        end
        set(H.spacefilepb,'Value',0);
        FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);
    end
    if FN.HRTFfiletype(1,1) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space_eq]));
    elseif FN.HRTFfiletype(1,1) == 2
        dir = 0;
        eval(['load -mat ' FN.space_path FN.space_eq]);
        TF1_space = TF1;    TF2_space = TF2;
        hrtfdirmat = dir;
        clear dir TF1 TF2
    else
        disp(['space3 HRTFfiletype incorrect'])
        return
    end
elseif XStimParams.ildalone_flag
    str1 = 'ILA';
    if isempty(FN.ildalone_eq)        % FN not yet picked
        [FN.ildalone_eq,FN.ILA_path] = uigetfile('*.*','Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.ildalonefile,'String',[FN.ILA_path FN.ildalone_eq]);
        end
        set(H.ildalonepb,'Value',0);
        FN.HRTFfiletype(2,1) = testHRTFfiletype(FN.ILA_path, FN.ildalone_eq);
    end
    if FN.HRTFfiletype(2,1) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone_eq]));
    elseif FN.HRTFfiletype(2,1) == 2
        dir = 0;
        eval(['load -mat ' FN.ILA_path FN.ildalone_eq]);
        TF1_ila = TF1;    TF2_ila = TF2;
        hrtfdirmat = dir; dir_ILA = dir;
        clear dir TF1 TF2
    else
        disp(['ildalone HRTFfiletype incorrect'])
        return
    end
elseif XStimParams.itdalone_flag
    str1 = 'ITA';
    if isempty(FN.itdalone_eq)        % FN not yet picked
        [FN.itdalone_eq,FN.ITA_path] = uigetfile('*.*','Select ITD-alone HRTF File');
        if(FN.ITA_path ~= 0)
            set(H.itdalonefile,'String',[FN.ITA_path FN.itdalone_eq]);
        end
        set(H.itdalonepb,'Value',0);
        FN.HRTFfiletype(3,1) = testHRTFfiletype(FN.ITA_path, FN.itdalone_eq);
    end
    if FN.HRTFfiletype(3,1) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.ITA_path FN.itdalone_eq]));
    elseif FN.HRTFfiletype(3,1) == 2
        dir = 0;
        eval(['load -mat ' FN.ITA_path FN.itdalone_eq]);
        TF1_ita = TF1;    TF2_ita = TF2;
        hrtfdirmat = dir;
        clear dir TF1 TF2
    else
        disp(['itdalone HRTFfiletype incorrect'])
        return
    end
    
    % also need FN.ildalone for this test
    if isempty(FN.ildalone_eq)        % FN not yet picked
        [FN.ildalone_eq,FN.ILA_path] = uigetfile('*.*','Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.ildalonefile,'String',[FN.ILA_path FN.ildalone_eq]);
        end
        set(H.ildalonepb,'Value',0);
        FN.HRTFfiletype(2,1) = testHRTFfiletype(FN.ILA_path, FN.ildalone_eq);
    end
    %%%%%%%
    if FN.HRTFfiletype(2,1) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone_eq]));
    elseif FN.HRTFfiletype(2,1) == 2
        dir = 0;
        eval(['load -mat ' FN.ILA_path FN.ildalone_eq]);
        TF1_ila = TF1;    TF2_ila = TF2;
        dir_ILA = dir;
        clear dir TF1 TF2
    else
        disp(['ildalone HRTFfiletype incorrect'])
        return
    end
end

% setup PD1
S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);

% reset DSPs
S232('PD1nstrms',1,2,2);		% 2 input and 2 output buffers
S232('PD1resetDSP',1,hex2dec('FFF'));
S232('dropall');
S232('PD1clrsched',1);

% make PD1 connections & load itd filters
if(XStimParams.ildalone_flag == 1) 					% ILDAlone Test
    % for left ear: IB[0] >> IREG[0] >> DSP[2] >> DSP[0] >> DAC[0]
    % for right ear:IB[1] >> IREG[1] >>  DSP[3] >> DSP[1] >> DAC[1]
    S232('PD1addsimp',1,S232('DSPout',0),S232('DAC',0));
    S232('PD1addsimp',1,S232('DSPout',1),S232('DAC',1));
    
    S232('PD1addsimp',1,S232('DSPout',2),S232('DSPin',0));
    S232('PD1addsimp',1,S232('DSPout',3),S232('DSPin',1));
    
    S232('PD1addsimp',1,S232('IREG',0),S232('DSPin',2));  
    S232('PD1addsimp',1,S232('IREG',1),S232('DSPin',3));  
    
    S232('PD1specIB',1,S232('IB',0),S232('IREG',0)); %IB to IREG
    S232('PD1specIB',1,S232('IB',1),S232('IREG',1));
    
    % calc ITD
    itdleft = 0; itdright = 0;
    if(XStimParams.curr_ITD < 0)
        itdright = abs(XStimParams.curr_ITD);
    elseif(XStimParams.curr_ITD > 0)
        itdleft = abs(XStimParams.curr_ITD);
    end
    
    %Load left ITDfilt to DSPs farthest from DACs   
    eval(['load ' FN.ITD_path 'itdfilt' num2str(itdleft)]);
    eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
    S232('pushf',itd_filt_left,length(itd_filt_left));
    S232('scale',.99);
    S232('PD1resetDSP',1,hex2dec('4'));
    S232('PreLoadRaw',1,S232('DSPid',2),'MONO','STACK','','',TDT.itd_scale,1.0,1);
    %Load right ITDfilt to DSPs farthest from DACs
    eval(['load ' FN.ITD_path 'itdfilt' num2str(itdright)]);
    eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
    S232('pushf',itd_filt_right,length(itd_filt_right));
    S232('scale',.99);
    S232('PD1resetDSP',1,hex2dec('8'));
    S232('PreLoadRaw',1,S232('DSPid',3),'MONO','STACK','','',TDT.itd_scale,1.0,1);
elseif(XStimParams.itdalone_flag == 1) 					% ITDAlone Test
    % for left ear: IB[0] >> IREG[0] >> DSP[2] >> DSP[0] >> DAC[0]
    % for right ear:IB[1] >> IREG[1] >>  DSP[3] >> DSP[1] >> DAC[1]
    S232('PD1addsimp',1,S232('DSPout',0),S232('DAC',0));
    S232('PD1addsimp',1,S232('DSPout',1),S232('DAC',1));
    
    S232('PD1addsimp',1,S232('DSPout',2),S232('DSPin',0));
    S232('PD1addsimp',1,S232('DSPout',3),S232('DSPin',1));
    
    S232('PD1addsimp',1,S232('IREG',0),S232('DSPin',2));  
    S232('PD1addsimp',1,S232('IREG',1),S232('DSPin',3));  
    
    S232('PD1specIB',1,S232('IB',0),S232('IREG',0)); %IB to IREG
    S232('PD1specIB',1,S232('IB',1),S232('IREG',1));
    
    %Load left MAGfilt to DSPs farthest from DACs
    if FN.HRTFfiletype(2,1) == 1
        dir_ILA = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone_eq]));
        ILAind = max(find(dir_ILA(1,:) == XStimParams.el & dir_ILA(2,:) == XStimParams.az));
        if isempty('ILAind')
            disp('Could not find ILA location in HRTF file');        return
        end
        eval(['ila_filt_left = mtlrch(''' FN.ILA_path FN.ildalone_eq ''', ' num2str(ILAind * 2-1) ');']);
        eval(['ila_filt_right = mtlrch(''' FN.ILA_path FN.ildalone_eq ''', ' num2str(ILAind * 2) ');']);
    else
        ILAind = max(find(dir_ILA(1,:) == XStimParams.el & dir_ILA(2,:) == XStimParams.az));
        if isempty('ILAind')
            disp('Could not find ILA location in HRTF file');        return
        end
        ila_filt_left = TF1_ila(ILAind);
        ila_filt_right = TF2_ila(ILAind);
    end
    S232('pushf',ila_filt_left,length(ila_filt_left));
    S232('scale',.99);
    S232('PD1resetDSP',1,hex2dec('4'));
    S232('PreLoadRaw',1,S232('DSPid',2),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
    %Load right MAGfilt to DSPs farthest from DACs
    S232('pushf',ila_filt_right,length(ila_filt_right));
    S232('scale',.99);
    S232('PD1resetDSP',1,hex2dec('8'));
    S232('PreLoadRaw',1,S232('DSPid',3),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
elseif(XStimParams.space_flag == 1) 			% Fully-cued Test
    % for left ear: IB[0] >> IREG[0] >> DSP[0] >>  DAC[0]
    S232('PD1addsimp',1,S232('DSPout',0),S232('DAC',0));
    S232('PD1addsimp',1,S232('IREG',0),S232('DSPin',0));  
    S232('PD1specIB',1,S232('IB',0),S232('IREG',0)); %IB to DSPin
    % for right ear:IB[1] >> IREG[0] >> DSP[1] >>  DAC[1]
    S232('PD1addsimp',1,S232('DSPout',1),S232('DAC',1));
    S232('PD1addsimp',1,S232('IREG',1),S232('DSPin',1));  
    S232('PD1specIB',1,S232('IB',1),S232('IREG',1));
end

%Specify DAMA buffers
clear BUF
BUF.L1				= 1;
BUF.R1				= 2;
BUF.L2				= 3;
BUF.R2				= 4;
BUF.playseq_L		= 5;
BUF.playseq_R		= 6;
BUF.playspec		= 7;

%Make PLAY sequence buffers
S232('dropall');
S232('allot16',BUF.playseq_L,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.L1);
S232('make',1,1);
S232('make',2,BUF.L2);
S232('make',3,1);
S232('make',4,0);
S232('qpop16',BUF.playseq_L);

S232('allot16',BUF.playseq_R,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.R1);
S232('make',1,1);
S232('make',2,BUF.R2);
S232('make',3,1);
S232('make',4,0);
S232('qpop16',BUF.playseq_R);

%Make play SPECification buffer
S232('allot16',BUF.playspec,11);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L);
S232('make',1,BUF.playseq_R);
S232('make',2,0);
S232('qpop16',BUF.playspec);

S232('allot16',BUF.L1, TDT.buf_pts);
S232('allot16',BUF.R1, TDT.buf_pts);
S232('allot16',BUF.L2, TDT.buf_pts);
S232('allot16',BUF.R2, TDT.buf_pts);

%Set MII parameters
mii_us_per_sample = 10; 				%microsecond per sample
mii_separation = 100; 					%only take events separated by 100 samples (i.e., 1 ms)

% set ABL
%S232('PA4mute',1);
S232('PA4mute',2);
S232('PA4atten',1,abs(XStimParams.curr_ABL)-20);
%S232('PA4atten',2,abs(XStimParams.curr_ABL)-20);

% get desired locations
if XStimParams.picker_flag				% use space picker
    XStimParams.locations = GUI.locations1';
else											% or read script
    locscriptfid = fopen([FN.script_path FN.script],'r');
    XStimParams.locations = fscanf(locscriptfid,'%i %i',[2 inf]);
    fclose(locscriptfid);
end

% write diary
if get(H.recorddata,'Value')
    tempstr = [str1 ' space3-test: ' num2str(length(XStimParams.locations)) ' locations;     atten: ' num2str(abs(XStimParams.curr_ABL))];
    update_diary
end

%Get the HRTF indices & HRTFs for desired locations
for locnum = 1:size(XStimParams.locations,2)
    if mod(locnum,100) == 1
        set(H.buildplay,'String','Loading HRTFs');
        if get(H.buildplay,'BackgroundColor') == [1 1 0]
            set(H.buildplay,'BackgroundColor','green');
        else
            set(H.buildplay,'BackgroundColor','yellow');
        end
    end
    
    locind(locnum) = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
    %read desired HRTFs
    if XStimParams.space_flag
        if FN.HRTFfiletype(1,1) == 1
            HRTF_L(:,locnum) = mtlrch([FN.space_path FN.space_eq],locind(locnum)*2-1);
            HRTF_R(:,locnum) = mtlrch([FN.space_path FN.space_eq],locind(locnum)*2);
        else
            HRTF_L(:,locnum) = TF1_space(locind(locnum),:)';
            HRTF_R(:,locnum) = TF2_space(locind(locnum),:)';
        end
    elseif XStimParams.ildalone_flag
        if FN.HRTFfiletype(1,1) == 1
            HRTF_L(:,locnum) = mtlrch([FN.ILA_path FN.ildalone_eq],locind(locnum)*2-1);
            HRTF_R(:,locnum) = mtlrch([FN.ILA_path FN.ildalone_eq],locind(locnum)*2);
        else
            HRTF_L(:,locnum) = TF1_ila(locind(locnum),:)';
            HRTF_R(:,locnum) = TF2_ila(locind(locnum),:)';
        end
    elseif XStimParams.itdalone_flag
        if FN.HRTFfiletype(1,1) == 1
            HRTF_L(:,locnum) = mtlrch([FN.ITA_path FN.itdalone_eq],locind(locnum)*2-1);
            HRTF_R(:,locnum) = mtlrch([FN.ITA_path FN.itdalone_eq],locind(locnum)*2);
        else
            HRTF_L(:,locnum) = TF1_ita(locind(locnum),:)';
            HRTF_R(:,locnum) = TF2_ita(locind(locnum),:)';
        end
    end
end      
set(H.buildplay,'BackgroundColor','yellow');
clear TF*

%%%%%% Make or get stimulus here
set(H.status,'String', 'Status: building stimulus');
%Make the stimuli
switch get(H.space3_stim_type,'Value')    
    case 1      % tone
        sound_left = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur,0);
    case 2         % GammaTones
        sound_left = MakeGammaTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur);
    case 6        %Narrowband Noise
        sound_left = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
    case 8        %Broadband Noise
        sound_left = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
    case 9        %Stimulus from file
        if isempty(FN.stim)                     % Stimulus from file
            [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select stimulus File');
        end
        fid = fopen([FN.stim_path FN.stim],'r');
        sound_left = fread(fid,inf,'float');
        fclose(fid);
    otherwise
        set(H.space3_stim_type,'Value',8);
        sound_left = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        disp('Unsupported stimulus type for space3 test; reset to BROADBAND');
end
% remove DC
sound_left = sound_left(:)' - mean(sound_left);
% scale to ACPower = 1;
sound_left = sound_left / mom(sound_left,2);
% R=L
sound_right = sound_left;

XStimParams.curr_stimdur = length(sound_left)* (1000/TDT.Fs);

% modulate
if get(H.space3_stim_type,'Value') ~= 9 & ~strcmp(XStimParams.mod_type,'None')
    Envelope = make_env(XStimParams.curr_stimdur, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
    sound_right = sound_right .* Envelope(:)';
    sound_left = sound_left .* Envelope(:)';
end

% remove DC
sound_left = sound_left - mean(sound_left);
sound_right = sound_right - mean(sound_right);

% Ramp the stimuli
if get(H.space3_stim_type,'Value') ~= 9
    ramp_time = 5; %ms
    [sound_left] = ramp_sound(sound_left,TDT.Fs,ramp_time);
    [sound_right] = ramp_sound(sound_right,TDT.Fs,ramp_time);
end

%Add silent period prior to stimulus to calculate spontaneous rate
silence_len = (DUR_silence * round(TDT.Fs/1000));
sound_left =  [zeros(1,silence_len) sound_left];
sound_right = [zeros(1,silence_len) sound_right];

%Add silent period after stimulus for HRTF filtering
sound_left =  [sound_left zeros(1,TDT.hrtffiltlen)];
sound_right = [sound_right zeros(1,TDT.hrtffiltlen)];

%Add silent period after stimulus for additional spike gathering
silence_len2 = (XStimParams.silence_trail * round(TDT.Fs/1000));
sound_left =  [sound_left zeros(1,silence_len2)];
sound_right = [sound_right zeros(1,silence_len2)];

% calc npts to play
npts_totalplay = silence_len + (XStimParams.curr_stimdur * round(TDT.Fs/1000))+ TDT.hrtffiltlen + silence_len2;

if(XStimParams.ildalone_flag == 1)     %Add silent period after stimulus for ITD filtering
    sound_left =  [sound_left zeros(1,TDT.itdfiltlen)];
    sound_right = [sound_right zeros(1,TDT.itdfiltlen)];
    npts_totalplay = npts_totalplay + TDT.itdfiltlen;
elseif(XStimParams.itdalone_flag == 1) %Add silent period after stimulus for ILA filtering
    sound_left =  [sound_left zeros(1,TDT.hrtffiltlen)];
    sound_right = [sound_right zeros(1,TDT.hrtffiltlen)];
    npts_totalplay = npts_totalplay + TDT.hrtffiltlen;
end

%scale stimuli
% 3/7/07
sound_left = sound_left * TDT.scaleFactor;
sound_right = sound_right * TDT.scaleFactor;

% %%%%%%% write sounds to disk
FNtempL = [FN.temp_stim_path 'L1'];
FNtempR = [FN.temp_stim_path 'R1'];
fid = fopen(FNtempL,'w');
fwrite(fid,sound_left,'short');
fclose(fid);
fid = fopen(FNtempR,'w');
fwrite(fid,sound_right,'short');
fclose(fid);
S232('PD1npts',1,npts_totalplay);	

pause(5);
set(H.buildplay,'String','Playing Stimuli');
set(H.buildplay,'BackgroundColor','yellow');
set(H.remreps,'String',num2str(XStimParams.numreps));

% increment testnum if this file exists
if(exist1('H.spacefig3') & get(H.recorddata,'Value'))
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%% play one blank sound first to work-around TDT glitch
if ~get(H.buildOnly,'value')
    set(H.status,'String', 'Status: playing blank trial');
S232('dropall');
% preload HRTFs to DSPs nearest to DACs
S232('PD1resetDSP',1, hex2dec('1'));
S232('pushf',HRTF_L(:,1)',255);
S232('PreLoadRaw',1,S232('DSPid',0),'MONO','STACK','','',TDT.hrtf_scale,TDT.hrtf_scale,1);
S232('PD1resetDSP',1, hex2dec('2'));
S232('pushf',HRTF_R(:,1)',255);
S232('PreLoadRaw',1,S232('DSPid',1),'MONO','STACK','','',TDT.hrtf_scale,TDT.hrtf_scale,1);
S232('seqplay',BUF.playspec);
%Set up MII
m100x( C_.INIT );
m110dx( C_.INIT );
m110dx( C_.CLOCK, mii_us_per_sample);
m110dx( C_.MODE, M_.PST );
%Fill first stimulus buffers
s_pt = 1;
S232('disk2dama16', BUF.L1,FNtempL,s_pt);
S232('disk2dama16', BUF.R1,FNtempR,s_pt);

%Fill second stimulus buffers
s_pt = (s_pt + TDT.buf_pts);
S232('disk2dama16', BUF.L1,FNtempL,s_pt);
S232('disk2dama16', BUF.R1,FNtempR,s_pt);


% set toggles
toggle1 = BUF.L1;
toggle2 = BUF.R1;

% set double buffering params
recordBUFIDs = [0 0 0 0];
decimationBUFID_factors = [0 0 0 0];
playBUFIDs = [BUF.L1 BUF.L2 BUF.R1 BUF.R2];
bufSTARTPOSITIONs = [0 0 (s_pt + TDT.buf_pts) (s_pt + TDT.buf_pts)];
bufLENGTHs = [0 0 TDT.buf_pts TDT.buf_pts];

% fill buffers
S232('dropall');
S232('dpush',TDT.buf_pts);
S232('tone',6000,1e6/TDT.Fs);
S232('scale',10000);
S232('qdup');
S232('qpop16',BUF.L1);
S232('qdup');
S232('qpop16',BUF.R1);

S232('tone',6000,1e6/TDT.Fs);
%S232('flat');
S232('scale',10000);
S232('qdup');
S232('qpop16',BUF.L2);
S232('qpop16',BUF.R2);
S232('dropall');

%S232('PD1arm',1);
%S232('pfireall');
%S232('PD1strig',1);

%Start clock
m110dx( C_.START);
%Send pulse: PD1 GO!
m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse

if 0
    %%%%% cycle through play buffers
    while S232('PD1status',1)
        s_pt = s_pt + TDT.buf_pts;
        %%%%% load next buffers
        if (s_pt <= npts_totalplay)
            if (toggle1 == BUF.L1)			% load second buffers
                S232('disk2dama16', BUF.L2,FNtempL,s_pt-1);
                S232('disk2dama16', BUF.R2,FNtempR,s_pt-1);
            else							% load first buffers
                S232('disk2dama16', BUF.L1,FNtempL,s_pt-1);
                S232('disk2dama16', BUF.R1,FNtempR,s_pt-1);
            end
        end
        %Wait for these buffers to finish playing
        while(S232('PD1status',1) &...
                S232('playseg',1) == toggle1 &...
                S232('playseg',2) == toggle2)
        end
        % toggle buffers
        if toggle1 == BUF.L1
            toggle1 = BUF.L2; toggle2 = BUF.R2;
        else
            toggle1 = BUF.L1; toggle2 = BUF.R1;
        end            
    end                                 %%%%% end of play buffer toggling loop
else    
    double_buffer('e:\kip\noisetc','e:\kip\noisetc', FNtempL,FNtempR, ...
        recordBUFIDs,decimationBUFID_factors,playBUFIDs, ...
        bufSTARTPOSITIONs,bufLENGTHs,npts_totalplay);    
end

%Stop the m110 and get spikes
m110dx( C_.STOP);
spikes = m110dx( C_.DATA, 1000); 			%Take 1000 spikes max

% set ABL
if(XStimParams.curr_ABL < -110) 
    XStimParams.curr_ABL = -110;
end
S232('PA4atten',1,abs(XStimParams.curr_ABL)-20);
%S232('PA4atten',2,abs(XStimParams.curr_ABL)-20);
% end of first blank sound


repnum = 1;
datamatrix = [];
numtrials = size(XStimParams.locations,2);
finalspikematrix = zeros(1,numtrials);

%%% start looping through reps 
while (exist1('H.spacefig3') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli for each rep
    randseq = randperm(numtrials);   
    trialnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    set(H.status,'String', ['Status: playing rep #' num2str(repnum)]);

    %%% start looping through trials (locations)
    tic
    while (exist1('H.spacefig3') & (trialnum <= numtrials))  
        S232('dropall');
        
        %Check for pause by user
        if pause_check  return; end
        
        % preload HRTFs to DSPs nearest to DACs
        %Load left HRTF
        S232('PD1resetDSP',1, hex2dec('1'));
        S232('pushf',HRTF_L(:,randseq(trialnum))',255);
        S232('PreLoadRaw',1,S232('DSPid',0),'MONO','STACK','','',TDT.hrtf_scale,TDT.hrtf_scale,1);
        %Load right HRTF
        S232('PD1resetDSP',1, hex2dec('2'));
        S232('pushf',HRTF_R(:,randseq(trialnum))',255);
        S232('PreLoadRaw',1,S232('DSPid',1),'MONO','STACK','','',TDT.hrtf_scale,TDT.hrtf_scale,1);
        S232('seqplay',BUF.playspec);
        
        %Set up MII
        m100x( C_.INIT );
        m110dx( C_.INIT );
        m110dx( C_.CLOCK, mii_us_per_sample);
        m110dx( C_.MODE, M_.PST );
        
        while toc < .1+XStimParams.test_ISI/1000 
            if pause_check    return; end
        end
        
        %Fill first stimulus buffers
        s_pt = 1;    
        S232('disk2dama16', BUF.L1,FNtempL,s_pt);
        S232('disk2dama16', BUF.R1,FNtempR,s_pt);
        
        %Fill second stimulus buffers
        s_pt = (s_pt + TDT.buf_pts);
        S232('disk2dama16', BUF.L1,FNtempL,s_pt-1);
        S232('disk2dama16', BUF.R1,FNtempR,s_pt-1);


        % set toggles
        toggle1 = BUF.L1;
        toggle2 = BUF.R1;
        
        S232('PD1arm',1);
        S232('pfireall');
        S232('PD1strig',1);
        
        if (trialnum <= numtrials)
            %Start clock
            m110dx( C_.START);
            %Send pulse: PD1 GO!
            m101x( C_.DATA,M_.BIT,M_.PULSE,0);  % Use port 0 for the pulse
        end
        tic
        
        if 0
            %%%%% cycle through play buffers
        while S232('PD1status',1)
            s_pt = s_pt + TDT.buf_pts;
            %%%%% load next buffers
            if (s_pt <= npts_totalplay)
                %if e_pt > npts_totalplay
                %    pts = npts_totalplay-s_pt+1;
                %else
                %    pts = buf_pts;
                %end
                if (toggle1 == BUF.L1)			% load second buffers
                    %S232('push16',sound_left(s_pt:s_pt+pts-1),pts);
                    %S232('qpop16',BUF.L2);
                    %S232('push16',sound_right(s_pt:s_pt+pts-1),pts);
                    %S232('qpop16',BUF.R2);
                    S232('disk2dama16', BUF.L2,FNtempL,s_pt-1);
                    S232('disk2dama16', BUF.R2,FNtempR,s_pt-1);
                else							% load first buffers
                    %S232('push16',sound_left(s_pt:s_pt+pts-1),pts);
                    %S232('qpop16',BUF.L1);
                    %S232('push16',sound_right(s_pt:s_pt+pts-1),pts);
                    %S232('qpop16',BUF.R1);
                    S232('disk2dama16', BUF.L1,FNtempL,s_pt-1);
                    S232('disk2dama16', BUF.R1,FNtempR,s_pt-1);
                end
            end       %%%%% finished loading buffers
            
            %Wait for these buffers to finish playing
            while(S232('PD1status',1) &...
                    S232('playseg',1) == toggle1 &...
                    S232('playseg',2) == toggle2)
            end
            % toggle buffers
            if toggle1 == BUF.L1
                toggle1 = BUF.L2;   toggle2 = BUF.R2;
            else
                toggle1 = BUF.L1;   toggle2 = BUF.R1;
            end            
        end                                 %%%%% end of play buffer toggling loo
    else
        double_buffer('e:\kip\noisetc','e:\kip\noisetc', FNtempL,FNtempR, ...
            recordBUFIDs,decimationBUFID_factors,playBUFIDs, ...
            bufSTARTPOSITIONs,bufLENGTHs,npts_totalplay);    
    end
    
    
        %Stop the m110 and get spikes
        m110dx( C_.STOP);
        spikes = m110dx( C_.DATA, XStimParams.curr_stimdur*2); 			%Take 2*XStimParams.curr_stimdur spikes max
        ind = find(spikes ~= 0); 	% Get clock events that are spikes
        spikes = spikes(ind);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; % Only take those events separated by >=1 ms
        
        % save spiketimes to each data array
        if(exist1('H.spacefig3') & ~isempty(spikes)) & (trialnum > 0)
            spikes = spikes(ind);
            spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
            EL_trial = [EL_trial;hrtfdirmat(1,locind(randseq(trialnum)))* ones(size(spikes))];
            AZ_trial = [AZ_trial;hrtfdirmat(2,locind(randseq(trialnum)))* ones(size(spikes))];
            repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
            Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
        end   
        
        %Check for pause by user
        if pause_check  return;     end
        
        if(trialnum > 0)
            finalspikematrix(randseq(trialnum)) = ...
                finalspikematrix(randseq(trialnum)) + length(spikes);
        end
        
        remtrials = numtrials - trialnum;
        set(H.remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end 
    %%%% end loop over trials (locations)
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.spacefig3') & ~exist1('H.finalspikeratefig'))
        H.finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','Space Test Spike Rate Plot',...
            'NumberTitle','off');
        H.spikeaxes = axes;
    end
    figure(H.finalspikeratefig)
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
    colorbar
    
    %Record Data
    if(exist1('H.spacefig3') & get(H.recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
        record_data3(XStimParams,datamatrix);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 
%%%%% end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.finalspikeratefig)
set(H.finalspikeratefig,'Name','Final Plot for Space Test');
plotdiam1(XStimParams.locations, interimspikerate);set(H.spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
        ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))],'FontSize',8);
colorbar
end     % if ~H.buildOnly

set(H.buildplay,'String','Build/Play status');
set(H.exitspace,'Visible','on');
set(H.resetspace,'Visible','on');

% increment test number
if(exist1('H.spacefig3') & get(H.recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
end
set(H.status,'String', 'Status: awaiting instructions');


%%%%%%%%%
function [flag] = pause_check

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Check for pause by user
flag = 0;
while (exist1('H.spacefig3') & get(H.pausespace,'Value'))
    pause(0);
    if(~exist1('H.spacefig3')) return; end         
    set(H.exitspace,'Visible','on');
    set(H.resetspace,'Visible','on');
    if(exist1('H.spacefig3') & get(H.resetspace,'Value') == 1)
        set(H.resetspace,'Value',0);
        set(H.pausespace,'Value',0);
        Reset_Space3;   flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_Space3;   flag = 1;
        return;
    end
end
if XStimParams.reset_flag
    flag = 1;
    XStimParams.reset_flag = 0;
end