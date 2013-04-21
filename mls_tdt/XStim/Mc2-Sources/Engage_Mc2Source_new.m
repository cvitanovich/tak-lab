function [] = Engage_Mc2Source()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI
global FigParams

%Engage_Mc2Source

%*******************************************************************************
%	The Mc2Source Test operation
%*******************************************************************************
% NEW version to use *.eq HRTFfiles
%       
% Jun08_09 mult by 1.732 to scale to mom2
% Jun08_09 changed default scaling to .5, so that total level will
% approximate one source

rand('state',sum(100*clock));

stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);

% check for normal or reversed S1/S2
if strcmp(XStimParams.test_type,'Mc2SourceReversed')
    flag_reverse = 1;
else
    flag_reverse = 0;
end

% Check for scaleFactor
scaleFactor = 60;
if ~XStimParams.HiDynamicRange
    BN = questdlg('Use HighDynamic Range?','Dynamic Range is set to low','High','Low (normal)','High');
    if ~isempty(findstr(BN,'Low'))    scaleFactor = TDT.scaleFactor;      end
end

set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

% calc stim duration
DUR = XStimParams.curr_stimdur;
nPts = DUR*round(TDT.Fs/1000);
nEpochmsecs = XStimParams.epoch_duration(1);       % duration of each epoch at one SPL
nEpochPts = round(nEpochmsecs * TDT.Fs/1000);
nEpochs = round(nPts / nEpochPts);
nPts = nEpochs * nEpochPts;
DUR = nPts / round(TDT.Fs/1000);
XStimParams.curr_stimdur = DUR;

% SPLs for low prob and for hi prob
SPLs_hi = 10+(XStimParams.focalSPL(1) :XStimParams.SPLinc(1): XStimParams.focalSPL(2));
SPLs_lo = setdiff([10:XStimParams.SPLinc(1):40], SPLs_hi);
if isempty(SPLs_lo)
    XStimParams.focalProb = 1;
    set(H.focalProb,'string','1');
end

% check if filt files assigned
if exist([FN.ephone_path FN.ephone2]) ~= 2
    ephonefilediagbox;
end

if XStimParams.space_flag
    while exist1([FN.space_path FN.space_eq]) ~=2
        [FN.space_eq,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select Fully-cued HRTF File *.eq preferred');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space_eq]);
        end
        FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);
    end
    disp('This is a FULLY-CUED Mc2source test')
else
    while exist1([FN.space_path FN.ablequal_eq]) ~=2
        [FN.ablequal_eq,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select ABLequal FC HRTF File  *.eq preferred');
        if(FN.space_path ~= 0)
            set(H.ABLequalfile,'String',[FN.space_path FN.ablequal_eq]);
        end
        FN.HRTFfiletype(7,1) = testHRTFfiletype(FN.space_path, FN.ABLequal_eq);
    end
    disp('This is an ABLequal FC Mc2source test')
end    

%Put parameters into XStimParams
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.reset_flag = 0;

if (get(H.stim_type,'Value') == 9) & isempty(FN.stim)
    [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select stimulus File');
end

%Specify DAMA buffers
clear BUF
BUF.snd1			= 1;
BUF.snd2			= 2;
BUF.playseq1		= 5;
BUF.playseq2		= 6;
BUF.playspec1		= 7;

%Make play sequence buffers
S232('allot16',BUF.playseq1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.snd1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq1);

S232('allot16',BUF.playseq2,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.snd2);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq2);

%Make play specification buffer
S232('allot16',BUF.playspec1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq1);
S232('make',1,BUF.playseq2);
S232('make',2,0);
S232('qpop16',BUF.playspec1);

%Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
silence_len = (XStimParams.silence_lead * round(TDT.Fs/1000));
%Add a piece of silence after stimulus 
silence_len2 = (XStimParams.silence_trail * round(TDT.Fs/1000));

%Make Stimulus buffers
nPts_totalplay = (silence_len2 + silence_len + nPts) + 2* TDT.hrtffiltlen;
S232('dpush',nPts_totalplay);
S232('value',0);
S232('allot16',BUF.snd1,nPts_totalplay);
S232('qpop16',BUF.snd1);
S232('dpush',nPts_totalplay);
S232('value',0);
S232('allot16',BUF.snd2,nPts_totalplay);
S232('qpop16',BUF.snd2);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,nPts_totalplay);

S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('FFF'));           % was '0xFFF'
S232('dropall');

%Make connections for left ear
sources(1) = S232('DSPout',0);          % snd1
sources(2) = S232('DSPout',2);          % snd2
scales(1) = XStimParams.factor(1);
scales(2) = XStimParams.factor(2);
S232('PD1addmult',1,sources,scales,2,S232('DAC',0));

%Make connections for right ear
sources(1) = S232('DSPout',1);      % snd1
sources(2) = S232('DSPout',3);      % snd2
scales(1) = XStimParams.factor(1);
scales(2) = XStimParams.factor(2);
S232('PD1addmult',1,sources,scales,2,S232('DAC',1));

S232('PD1addsimp',1,S232('IREG',0),S232('DSPin',0));
S232('PD1addsimp',1,S232('IREG',0),S232('DSPin',1));
S232('PD1specIB',1,S232('IB',0),S232('IREG',0));

S232('PD1addsimp',1,S232('IREG',1),S232('DSPin',2));
S232('PD1addsimp',1,S232('IREG',1),S232('DSPin',3));
S232('PD1specIB',1,S232('IB',1),S232('IREG',1));

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4mute',1);	
S232('PA4mute',2);	

ISI = XStimParams.test_ISI;

%Get HRTF index
XStimParams.locations = GUI.locations1';
if size(XStimParams.locations,2) ~= 1
    set(H.pickerfig,'Color', [.1 .2 .8]);
    set(H.picker_error,'visible','on');
    disp('Mc2Source requires one and only one picked location')
    pause;
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end

%%%% load HRTFs
dir = 0;
if XStimParams.space_flag
    tempFN = FN.space_eq;
    filetype = FN.HRTFfiletype(1,1);
else
    tempFN = FN.ablequal_eq;
    filetype = FN.HRTFfiletype(7,1);
end
if filetype == 2
    eval(['load -mat ' FN.space_path tempFN]);
elseif filetype == 1
    dir = sph2dbl(mtlrdir([FN.space_path tempFN]));
end
XStimParams.HRTF_FN = tempFN;

% find locations for 2 sources NOTE:: OFFSETS here are for loc2 vs loc1
% (loc2 = loc1 + offset)
% (not, as in Two Source, either source vs mid-array)
offset_el = XStimParams.offset_el(2);
offset_az = XStimParams.offset_az(2);

loc1 = max(find(dir(1,:) == XStimParams.locations(1) & dir(2,:) == XStimParams.locations(2)));
loc2 = max(find(dir(1,:) == (XStimParams.locations(1) + offset_el) & dir(2,:) == (XStimParams.locations(2)) + offset_az));

str1 = ['  ' num2str(XStimParams.locations(1)) ' el, ' num2str(XStimParams.locations(2)) ' az, vs: ' ...
        num2str(XStimParams.locations(1) + offset_el) ' el, ' num2str(XStimParams.locations(2) + offset_az) ' az'];

if isempty(loc1)  || isempty(loc2)
    error('offset not possible, repick')
    return;
end

%load DSPs for sound#1 HRTF filtering
if FN.HRTFfiletype(1,1) == 1
    hrtf_left = mtlrch([FN.space_path XStimParams.HRTF_FN], 2*loc1 -1);
    hrtf_right = mtlrch([FN.space_path XStimParams.HRTF_FN], 2*loc1);
else
    hrtf_left = TF1(loc1,:);
    hrtf_right = TF2(loc1,:);
end
%Load left   
S232('pushf',hrtf_left,length(hrtf_left));
S232('PreLoadRaw',1,S232('DSPid',0),'MONO','STACK','','',TDT.ephonescale,1.0,1);
%Load right
S232('pushf',hrtf_right,length(hrtf_right));
S232('PreLoadRaw',1,S232('DSPid',1),'MONO','STACK','','',TDT.ephonescale,1.0,1);

%load DSPs for sound#2 HRTF filtering
if FN.HRTFfiletype(1,1) == 1
    hrtf_left = mtlrch([FN.space_path XStimParams.HRTF_FN], 2*loc2 -1);
    hrtf_right = mtlrch([FN.space_path XStimParams.HRTF_FN], 2*loc2);
else
    hrtf_left = TF1(loc2,:);
    hrtf_right = TF2(loc2,:);
end
%Load left   
S232('pushf',hrtf_left,length(hrtf_left));
S232('PreLoadRaw',1,S232('DSPid',2),'MONO','STACK','','',TDT.ephonescale,1.0,1);
%Load right
S232('pushf',hrtf_right,length(hrtf_right));
S232('PreLoadRaw',1,S232('DSPid',3),'MONO','STACK','','',TDT.ephonescale,1.0,1);

% write diary
if get(H.recorddata,'Value')
    str3 = '';
    if flag_reverse
        str0 = '     Mc2SrcRev  ';
        str2 = ['   SPL2:' num2str(10+XStimParams.source1_SPL) ...
                '  SPL1:' num2str(10+XStimParams.focalSPL(1)) '-' num2str(10+XStimParams.focalSPL(2))];
        
        if loc1 == loc2     % SPL1 and SPL2 from same loc
           str1 = ['  ' num2str(XStimParams.locations(1)) ' el, ' num2str(XStimParams.locations(2)) ' az'];
           if scales       % IG
                str3 = ' IG    ';
            else
                str3 = ' noG    ';
                str2 = ['  SPL1:' num2str(10+XStimParams.focalSPL(1)) '-' num2str(10+XStimParams.focalSPL(2))];
            end
        else
            str3 = ' DG    ';
        end
    else
        str0 = '     Mc2Src';
        str2 = ['   SPL1(' num2str(XStimParams.epoch_duration(1)) '):' num2str(10+XStimParams.source1_SPL) ...
                '  SPL2(' num2str(XStimParams.epoch_duration(2)) '):' num2str(10+XStimParams.focalSPL(1)) '-' num2str(10+XStimParams.focalSPL(2))];
    end
    tempstr = [str0 str3 str1 str2 ...
            ' p = ' num2str(XStimParams.focalProb) '   atten: ' num2str(abs(XStimParams.curr_ABL)) 'dB'];
    update_diary
end
clear TF1 TF2 dir hrtf_* loc1 loc2

%%%%%%%% play out a quiet sound while stims are being made
tic
S232('PA4mute',1);	
S232('PA4mute',2);	

S232('seqplay',BUF.playspec1);
S232('PD1arm',1);

%Send trigger
%Set up MII
m100x( C_.INIT );
m110dx( C_.INIT );
m110dx( C_.CLOCK, mii_us_per_sample);
m110dx( C_.MODE, M_.PST );

%Start clock & Send pulse: PD1 GO!
tic
m110dx( C_.START);
m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse

%%%%%%%%%%%%%%%%%%%% make the stimuli we'll use
set(H.status,'BackgroundColor','blue');
set(H.status,'String','Status: Building Stimuli');
set(H.remreps,'String','1');
finalspikematrix = 0;

% make SPL sequences for each source
% source2 will jitter around in SPL, but fall to 0 SPL whenever source1 comes on
if isempty(SPLs_lo)
    XStimParams.focalProb = 1;
end

if XStimParams.repeatmode
    dur(1) = XStimParams.onEpochs + XStimParams.offEpochs+1;        % total epochs of each rep
    n(1) = XStimParams.onEpochs/(XStimParams.source1_ISI);        % #probes during adapter ON
    n(2) = XStimParams.offEpochs/(XStimParams.source1_ISI);       % #probes during adapter OFF
else
    dur(1) = nEpochs;        % total epochs of each rep
    n(1) = round(nEpochs/(XStimParams.source1_ISI));        % #probes during adapter ON    
end
temp_ones = ones(nEpochPts,1);

% to hold AC sequences for each location and each rep
ACseq1 = cell(XStimParams.numreps,1);
ACseq2 = ACseq1;

% make SPLseq
randEpochs = rand(nEpochs+10,1);        
randEpochs2 = round(rand(nEpochs+10,1)*length(SPLs_hi)+.49999999999);  
randEpochs3 = round(rand(nEpochs+10,1)*length(SPLs_lo)+.49999999999);  

seq1 = zeros(nEpochs,1);                             % SPL for each epoch
seq2 = zeros(nEpochs,1);                             % SPL for each epoch
seqSPL1 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt
seqSPL2 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt

rampEnv = ones(nEpochPts,1);
if XStimParams.ramp
    rampEnv(1:150) = linspace(0,1,150);
    rampEnv(end-149:end) = linspace(1,0,150);
end

for iepoch = 1:dur(1):nEpochs
    % probes
    ind1 = [];
    if XStimParams.probeFlag(1)
        ind1 = [iepoch:XStimParams.source1_ISI:iepoch+XStimParams.source1_ISI*(n(1))];
    end
    if XStimParams.repeatmode & XStimParams.probeFlag(2)
        ind1 = [ind1 (iepoch+XStimParams.source1_ISI*(1+n(1))-1):XStimParams.source1_ISI:iepoch+XStimParams.source1_ISI*(n(1)+n(2))];
    end
    
    for jepoch = 1:length(ind1)
        ind3 = ind1(jepoch);
        if flag_reverse     % in reverse mode, probes vary in level, adapter is constant
            if randEpochs(jepoch+1) <= XStimParams.focalProb
                seq1(ind3) = SPLs_hi(randEpochs2(jepoch+1));
            else
                seq1(ind3) = SPLs_lo(randEpochs3(jepoch+1));
            end
        else
            seq1(ind3) = 10+XStimParams.source1_SPL;
        end
        ind0 = (ind3-1)*nEpochPts+1:ind3*nEpochPts;
        
        if flag_reverse 
            seqSPL1(ind0) = seq1(ind1(jepoch)) * rampEnv;
        else
            seqSPL1(ind0) = seq1(ind1(jepoch)) * ones(size(ind0));
        end
    end
    
    % adapter on during onEpochs when probe is not on
    if XStimParams.repeatmode
        ind2 = setdiff([iepoch:iepoch+XStimParams.onEpochs-1],ind1);
    else
        ind2 = setdiff([1:dur(1)],ind1);
    end
    for jepoch = 1:length(ind2) 
        ind3 = ind2(jepoch);
        if flag_reverse     % in reverse mode, probes vary in level, adapter is constant
            seq2(ind3) =  10+XStimParams.source1_SPL;
        else
            if randEpochs(ind3) <= XStimParams.focalProb
                seq2(ind3) = SPLs_hi(randEpochs2(ind3));
            else
                seq2(ind3) = SPLs_lo(randEpochs3(ind3));
            end
        end
        ind0 = (ind3-1)*nEpochPts+1:ind3*nEpochPts;
        if flag_reverse 
            seqSPL2(ind0) = seq2(ind2(jepoch)) * ones(size(ind0));
        else
            seqSPL2(ind0) = seq2(ind2(jepoch)) * rampEnv;
        end
    end
end   

clear ind0 ind1 ind2 ind3


% dB to linear
ACseq1{1} = seq1(:)';             % not updated - same for every rep
ACseq2{1} = seq2(:)';
seqSPL1 = 1.732 * 10.^(seqSPL1/20)-1;       % Jun08_09 mult by 1.732 to scale to mom2
seqSPL2 = 1.732 * 10.^(seqSPL2/20)-1;       % Jun08_09 mult by 1.732 to scale to mom2
if length(seqSPL1) ~= nPts
    disp(['Euston, we have a problem ...']);
end    

%Check for pause by user
if pause_check  return; end

%%%%%%%%%%%%%%%%%%%%%% make stim1 and load to buffer
param0 = get(H.stim_type,'Value');
if param0 ~= 9
    stim = get_stim(param0,XStimParams.freq(1));
else
    stim = get_stim(param0,FN.stim_path,FN.stim);
end

% remove any DCoffset
stim = stim - mom(stim,1);

% scale to AC power == 1
stim = stim/mom(stim,2);

% modulate stim1
if ~strcmp(XStimParams.mod_type,'None')
    Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
    stim = stim .* Envelope(:)';
end

%Ramp stim1
if get(H.stim_type,'Value') ~= 9
    ramp_time = 5; %ms
    [stim] = ramp_sound(stim,TDT.Fs,ramp_time);
end

% remove any DCoffset and weight
stim = (stim - mom(stim,1)) * XStimParams.factor(1) * scaleFactor;

%%%%% multiply by the epoch-envelopes for loc1
%lag = silence_len+100;          % +100 is an approximation of the HRTF-delay
%stim(lag:lag+nPts-1) = stim1 .* seqSPL1;

S232('push16',[zeros(1,silence_len) stim.* seqSPL1 zeros(1,silence_len2+ 2* TDT.hrtffiltlen)],nPts_totalplay);
S232('qpop16',BUF.snd1);

%%%% make Stim 2 and load to buffer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param0 = get(H.stim_type2,'Value');
if param0 ~= 9
    stim = get_stim(param0,XStimParams.freq(2));
else
    stim = get_stim(param0,FN.stim_path2,FN.stim2);
end
% remove any DCoffset
stim = stim - mom(stim,1);
% scale to ACpower == 1
stim = stim/mom(stim,2);

% modulate stim2
if ~strcmp(XStimParams.mod_type2,'None')
    Envelope = make_env(DUR, XStimParams.mod_type2, XStimParams.mod_depth(2), XStimParams.mod_freq(2), XStimParams.mod_phase(2));
    stim = stim .* Envelope(:)';
end

%Ramp the stimuli
if get(H.stim_type2,'Value') ~= 9
    ramp_time = 5; %ms
    stim = ramp_sound(stim,TDT.Fs,ramp_time);
end

% remove any DCoffset and weight
stim = (stim - mom(stim,1)) * XStimParams.factor(2) * scaleFactor;

%%%%% multiply by the epoch-envelopes for loc2
S232('push16',[zeros(1,silence_len) stim.* seqSPL2 zeros(1,silence_len2+ 255*2)],nPts_totalplay);
S232('qpop16',BUF.snd2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% finished making sounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);
pause(.1)        % in case there is a click on atten reset

%%%%%%%%%%%%% WAIT HERE FOR TOC


% Begin playing sounds   
set(H.status,'String','Status: Playing Stimuli');
set(H.status,'BackgroundColor','green');
set(H.remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.Mc2Sourcefig') & get(H.recorddata,'Value'))
    update_dataFN;
end
firstloop=1;
remtrials = XStimParams.numreps;
%loop for reps
while (exist1('H.Mc2Sourcefig') & (repnum <= XStimParams.numreps))
    trialnum = 1-firstloop;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    %Check for pause by user
    if pause_check  return; end
    
    S232('seqplay',BUF.playspec1);
    S232('PD1arm',1);
    %Set up MII & send trigger
    m100x( C_.INIT );
    m110dx( C_.INIT );
    m110dx( C_.CLOCK, mii_us_per_sample);
    m110dx( C_.MODE, M_.PST );
    
    while toc < (XStimParams.curr_stimdur + ISI)/1000     end
    tic
    m110dx( C_.START);
    m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
    
    if repnum < XStimParams.numreps     % make next stim
        % make new SPLseq
        randEpochs = rand(nEpochs+10,1);                                        % to choose between hiSPLs or loSPLs  
        randEpochs2 = round(rand(nEpochs+10,1)*length(SPLs_hi)+.49999999999);   % to choose which hi SPL 
        randEpochs3 = round(rand(nEpochs+10,1)*length(SPLs_lo)+.49999999999);   % to choose which lo SPL
        
        seq1 = zeros(nEpochs,1);                             % SPL for each epoch
        seq2 = zeros(nEpochs,1);                             % SPL for each epoch
        seqSPL1 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt
        seqSPL2 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt
        
        for iepoch = 1:dur(1):nEpochs
            % probes
            ind1 = [];
            if XStimParams.probeFlag(1)
                ind1 = [iepoch:XStimParams.source1_ISI:iepoch+XStimParams.source1_ISI*(n(1))];
            end
            if XStimParams.repeatmode
                if XStimParams.probeFlag(2)
                    ind1 = [ind1 (iepoch+XStimParams.source1_ISI*(1+n(1))-1):XStimParams.source1_ISI:iepoch+XStimParams.source1_ISI*(n(1)+n(2))];
                end
            end
            
            for jepoch = 1:length(ind1)
                ind3 = ind1(jepoch);
                if flag_reverse     % in reverse mode, probes vary in level, adapter is constant
                    if randEpochs(jepoch+1) <= XStimParams.focalProb
                        seq1(ind3) = SPLs_hi(randEpochs2(jepoch+1));
                    else
                        seq1(ind3) = SPLs_lo(randEpochs3(jepoch+1));
                    end
                else
                    seq1(ind3) = 10+XStimParams.source1_SPL;
                end
                ind0 = (ind3-1)*nEpochPts+1:ind3*nEpochPts;
                
                if flag_reverse 
                    seqSPL1(ind0) = seq1(ind1(jepoch)) * rampEnv;
                else
                    seqSPL1(ind0) = seq1(ind1(jepoch)) * ones(size(ind0));
                end
            end
            
            % adapter on during onEpochs when probe is not on
            if XStimParams.repeatmode
                ind2 = setdiff([iepoch:iepoch+XStimParams.onEpochs-1],ind1);
            else
                ind2 = setdiff([1:dur(1)],ind1);
            end
            for jepoch = 1:length(ind2) 
                ind3 = ind2(jepoch);
                if flag_reverse     % in reverse mode, probes vary in level, adapter is constant
                    seq2(ind3) =  10+XStimParams.source1_SPL;
                else
                    if randEpochs(ind3) <= XStimParams.focalProb
                        seq2(ind3) = SPLs_hi(randEpochs2(ind3));
                    else
                        seq2(ind3) = SPLs_lo(randEpochs3(ind3));
                    end
                end
                ind0 = (ind3-1)*nEpochPts+1:ind3*nEpochPts;
                if flag_reverse 
                    seqSPL2(ind0) = seq2(ind2(jepoch)) * ones(size(ind0));
                else
                    seqSPL2(ind0) = seq2(ind2(jepoch)) * rampEnv;
                end
            end
        end   
        clear ind0 ind1 ind2 ind3
        
        % dB to linear
        ACseq1{repnum+1} = seq1(:)';             % not updated - same for every rep
        seqSPL1 = 1.732 * 10.^(seqSPL1/20)-1;       % Jun08_09 mult by 1.732 to scale to mom2
        ACseq2{repnum+1} = seq2(:)';
        seqSPL2 = 1.732 * 10.^(seqSPL2/20)-1;       % Jun08_09 mult by 1.732 to scale to mom2
        if length(seqSPL1) ~= nPts
            disp(['Euston, we have a problem ...']);
        end    
        
        %Check for pause by user
        if pause_check  return; end
        
        %%% make stim1
        if param0 ~= 9
            stim = get_stim(param0,XStimParams.freq(1));
        else
            stim = get_stim(param0,FN.stim_path,FN.stim);
        end
        % remove any DCoffset
        stim = stim - mom(stim,1);
        % scale to AC power == 1
        stim = stim/mom(stim,2);
        % modulate stim1
        if ~strcmp(XStimParams.mod_type,'None')
            Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
            stim = stim .* Envelope(:)';
        end
        %Ramp stim1
        if get(H.stim_type,'Value') ~= 9
            ramp_time = 5; %ms
            [stim] = ramp_sound(stim,TDT.Fs,ramp_time);
        end
        % remove any DCoffset and weight
        stim = (stim - mom(stim,1)) * XStimParams.factor(1) * scaleFactor;
        
        %%%%% multiply by the epoch-envelopes for loc1 and load to buffer
        S232('push16',[zeros(1,silence_len) stim.* seqSPL1 zeros(1,silence_len2+ 255*2)],nPts_totalplay);
        S232('qpop16',BUF.snd1);
        
        %%%% make Stim 2 and load to buffer
        param0 = get(H.stim_type2,'Value');
        if param0 ~= 9
            stim = get_stim(param0,XStimParams.freq(2));
        else
            stim = get_stim(param0,FN.stim_path2,FN.stim2);
        end
        % remove any DCoffset
        stim = stim - mom(stim,1);
        % scale to ACpower == 1
        stim = stim/mom(stim,2);
        % modulate stim2
        if ~strcmp(XStimParams.mod_type2,'None')
            Envelope = make_env(DUR, XStimParams.mod_type2, XStimParams.mod_depth(2), XStimParams.mod_freq(2), XStimParams.mod_phase(2));
            stim = stim .* Envelope(:)';
        end
        %Ramp the stimuli
        if get(H.stim_type2,'Value') ~= 9
            ramp_time = 5; %ms
            stim = ramp_sound(stim,TDT.Fs,ramp_time);
        end
        % remove any DCoffset and weight
        stim = (stim - mom(stim,1)) * XStimParams.factor(2) * scaleFactor;
        
        %%%%% multiply by the epoch-envelopes for loc2
        S232('push16',[zeros(1,silence_len) stim.* seqSPL2 zeros(1,silence_len2+ 255*2)],nPts_totalplay);
        S232('qpop16',BUF.snd2);

        if pause_check  return; end
        
    end             % end of making next stim
    
    %Wait for TOC
    while toc <= (XStimParams.curr_stimdur+2000)/1000     end
    S232('PD1stop',1);
    
    %Stop the m110 and get spikes
    m110dx( C_.STOP);
    spikes = m110dx( C_.DATA, DUR/2); 			% Take DUR/2 spikes max
    ind = find(spikes ~= 0); 						% Get clock events that are spikes
    spikes = spikes(ind);
    ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
    if ~isempty(spikes) 
        spikes = spikes(ind);
        spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
        EL_trial = [EL_trial;XStimParams.locations(1)* ones(size(spikes))];
        AZ_trial = [AZ_trial;XStimParams.locations(2)* ones(size(spikes))];
        repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
        Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
    end
    
    %Check for pause by user
    if pause_check  return; end
    
    %Plot Spike Rate Data
    finalspikematrix = finalspikematrix + length(spikes);
    interimspikerate = finalspikematrix/repnum;
    if ~exist1('H.finalspikeratefig')
        H.finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','Mc2Source Test Spike Rate Plot',...
            'NumberTitle','off');
        figure(H.finalspikeratefig);
        H.spikeaxes = axes;
    end
    figure(H.finalspikeratefig)
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
    colorbar
    
    % calc and plot RLF data
    if ~exist1('H.Mc2SourceRLFfig')
        H.Mc2SourceRLFfig = figure;
    end
    figure(H.Mc2SourceRLFfig)
    xlabel('SPL (dB)')
    ylabel('spike rate (sp/s)')
    title([FN.data ' response to probe as RLF, rep#' num2str(repnum)])
    if repnum == 1
        SPLbins = 0:40;
        RLFall = cell(size(SPLbins'));
        n1 = zeros(size(SPLbins'));
    end
    [RLFall, n1, SPLx, RLFm] = plot_Mc2srcA(spikes_trial, ACseq1{repnum}, RLFall, n1, repnum);            % plots the RLF
    
    %Record Data
    if get(H.recorddata,'Value')
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
        record_data3(XStimParams,datamatrix, [ACseq1 ACseq2]);
    end
    
    set(H.remreps,'String',num2str(XStimParams.numreps - repnum));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.finalspikeratefig)
set(H.finalspikeratefig,'Name','Final Plot for Mc2Source Test');
plotdiam1(XStimParams.locations, interimspikerate);
set(H.spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
        ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))], 'FontSize',8);
colorbar
S232('PA4mute',1);	
S232('PA4mute',2);	

% plot RLF on persistant Figure
if isempty(FigParams.handle)
    FigParams.handle = figure; hold on;
end
prompt{1} = 'Figure Handle for plotting RLF';
prompt{2} = 'marker symbol and color';
default{1} = num2str(FigParams.handle);
default{2} = 'bo-';
answer = inputdlg(prompt,'Plotting params for RLF',1,default);
FigParams.handle = eval(answer{1});
FigParams.marker = answer{2};
figure(FigParams.handle)
[temp ind0] = intersect(SPLbins,10+(XStimParams.focalSPL(1):XStimParams.focalSPL(2)));     % index into focal SPLs
[temp ind1] = setdiff(SPLbins(11:end),temp);            % index to non-focal SPLs
plot(SPLx(ind0) , RLFm(ind0),FigParams.marker(2:end),'color',FigParams.marker(1), 'MarkerFaceColor', FigParams.marker(1));
plot(SPLx(ind1+10) , RLFm(ind1+9),FigParams.marker(1:2));

set(H.status,'String','Status: Results');
set(H.exitMc2Source,'Visible','on');
set(H.resetMc2Source,'Visible','on');

% increment test number
if(exist1('H.Mc2Sourcefig') & get(H.recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stim] = get_stim(param0,param1, param2)
% param0: stim_type (value)
% param1: XStimParams.freq or FN.stim_path
% param2: FN.stim

global H
global XStimParams
global TDT

switch param0      
    case 1    %tone at specified frequency 
        stim = MakeTone(TDT.Fs,param1,XStimParams.curr_stimdur);
    case 2    %GammaTones
        stim = MakeGammaTone(TDT.Fs,param1,XStimParams.curr_stimdur);
    case 8    %Broadband Noise
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
    case 9    %Stimulus from file
        if(~exist1('stim_from_file'))
            fid = fopen([param1 param2],'r');
            stim = fread(fid,inf,'float');
            fclose(fid);
        end
        if(size(stim,1) > 1)
            stim = stim';
        end
     case 10                %Broadband Noise then flattened
         stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
         %stim = flatten(stim);
         env = abs(hilbert(stim));
         Ienv = max(env)./ env;
         stim = stim .* Ienv;
   otherwise
        set(H.stim_type,'Value',8);
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        disp('Stimulus type not supported for Mc2Source tests.  Reset to BROADBAND');
        return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Envelope] = make_env(DUR,mod_type,param1,param2,param3)

% param1 used as follows:
%    'tone' or 'LP noise' or 'File': mod_depth
% param2 used as follows:
%    'tone or 'LP noise' : mod_freq
% param3 used as follows
%    'tone' : addition to modulation starting phase (0 - pi)

global H
global TDT
global FN

if nargin < 3   param3 = 0; end

Npts = DUR*(round(TDT.Fs/1000));
Inc = 1/TDT.Fs;
switch mod_type
    case 'Tone'
        T = 0:Inc:(DUR/1000 - Inc);
        Tone = (param1 / 2)* sin(2 * pi * param2 .* T + (.75 * 2 * pi + param3));
        Envelope = Tone + (1-param1/2);
    case 'LP Noise'
        LP_noise = m_noi(5, param2, Npts/(TDT.Fs), TDT.Fs/2);
        LP_noise = (param1 / 2)* (LP_noise / max1(LP_noise));
        Envelope = LP_noise + (1-param1/2);
    case 'File'			
        if(~exist1('mod_from_file'))
            fid = fopen(FN.mod,'r');
            mod_from_file = fread(fid,inf,'float');
            fclose(fid);
            while length(mod_from_file) ~= Npts
                [FN.mod,FN.mod_path] = uigetfile('*.*','Select Envelope File');
                if(FN.mod_path ~= 0)
                    set(H.modfile,'String',[FN.mod_path FN.mod]);
                end
                fid = fopen(FN.mod,'r');
                mod_from_file = fread(fid,inf,'float');
                fclose(fid);
            end
            mod_from_file = mod_from_file - mean(mod_from_file);
            mod_from_file = (param1 / 2)* (mod_from_file / max1(mod_from_file));
        end
        Envelope = mod_from_file + (1-param1/2);
    otherwise
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flag] = pause_check

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

flag = 0;
%Check for pause by user
while (exist1('H.Mc2Sourcefig') & get(H.pauseMc2Source,'Value'))
    pause(0);
    if(~exist1('H.Mc2Sourcefig')) return; end         
    set(H.exitMc2Source,'Visible','on');
    set(H.resetMc2Source,'Visible','on');
    if(exist1('H.Mc2Sourcefig') & get(H.resetMc2Source,'Value') == 1)
        set(H.resetMc2Source,'Value',0);
        set(H.pauseMc2Source,'Value',0);
        Reset_Mc2Source;   flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_Mc2Source;   flag=1;
        return;
    end
end

if XStimParams.reset_flag
    flag = 1;
    XStimParams.reset_flag = 0;
end