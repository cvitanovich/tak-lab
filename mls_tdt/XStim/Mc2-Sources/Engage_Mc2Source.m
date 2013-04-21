function [] = Engage_Mc2Source()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_Mc2Source

%*******************************************************************************
%	The Mc2Source Test operation
%*******************************************************************************
% NOTE: (Mar 4, 2003) problem with writing stim files discovered:
%   locind has multiple rows and this was not taken into account when writing files
%   so actual file written (and played out) was ceil(locind(randseq(trialnum))/2)
%   
% NOTE: (Jun 1, 2003) problem with make_env that always used mod_freq and mod_depth
%   from stimulus number1 (so all stims had same mod depth and freq
%   changed to use input param1,param2 to carry the correct depth and freq
%   also: eliminated 2nd buffercycle and associated buffers...
%
% NOTE: (July 29, 2004) added in ABL-alone and removed filecache, itd2,3 and az2,3 and
% el2,3 assuming these will always be the same as #1
% removed location script ability (forcing spacepicker)

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
SPLs_hi = 10+(XStimParams.focalSPL(1) :XStimParams.SPLinc: XStimParams.focalSPL(2));
SPLs_lo = setdiff([10:XStimParams.SPLinc:40], SPLs_hi);
if isempty(SPLs_lo)
    XStimParams.focalProb = 1;
    set(H.focalProb,'string','1');
end


% check if filt files assigned
if exist([FN.ephone_path FN.ephone2]) ~= 2
    ephonefilediagbox;
end

if XStimParams.space_flag
    while exist1([FN.space_path FN.space_std]) ~=2
        [FN.space_std,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select Fully-cued HRTF File');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space_std]);
        end
        FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
    end
    disp('This is a FULLY-CUED Mc2source test')
else
    while exist1([FN.space_path FN.ablequal_std]) ~=2
        [FN.ablequal_std,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select ABLequal FC HRTF File');
        if(FN.space_path ~= 0)
            set(H.ABLequalfile,'String',[FN.space_path FN.ablequal_std]);
        end
        FN.HRTFfiletype(7,2) = testHRTFfiletype(FN.space_path, FN.ABLequal_std);
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
BUF.L1				= 1;
BUF.R1				= 2;
BUF.playseq_L		= 5;
BUF.playseq_R		= 6;
BUF.playspec1		= 7;

%Make play sequence buffers
S232('allot16',BUF.playseq_L,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.L1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_L);

S232('allot16',BUF.playseq_R,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.R1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_R);

%Make play specification buffer
S232('allot16',BUF.playspec1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L);
S232('make',1,BUF.playseq_R);
S232('make',2,0);
S232('qpop16',BUF.playspec1);

%Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
silence_len = (XStimParams.silence_lead * round(TDT.Fs/1000));
%Add a piece of silence after stimulus 
silence_len2 = (XStimParams.silence_trail * round(TDT.Fs/1000));

%Make Stimulus buffers
nPts_totalplay = (silence_len2 + silence_len + nPts) + TDT.ephonefiltlen + TDT.hrtffiltlen;
S232('dpush',nPts_totalplay);
S232('value',0);
S232('allot16',BUF.L1,nPts_totalplay);
S232('qpop16',BUF.L1);
S232('dpush',nPts_totalplay);
S232('value',0);
S232('allot16',BUF.R1,nPts_totalplay);
S232('qpop16',BUF.R1);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,nPts_totalplay);

%Get Earphone filters
if FN.HRTFfiletype(6) == 1
    ephonefname = [FN.ephone_path FN.ephone2];
    ephonefilt_left  = (mtlrch(ephonefname,1))';
    ephonefilt_right = (mtlrch(ephonefname,2))';
elseif FN.HRTFfiletype(6) == 2
    eval(['load ' FN.ephone_path FN.ephone2]);
    ephonefilt_left  = TF1;
    ephonefilt_right = TF2;
else
    disp('ephone filters of incorrect type');
end

%Load Earphone filters
dspid_left = 0; dspid_right = 1;
S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('FFF'));           % was '0xFFF'
S232('dropall');
%Make connections for left ear
S232('PD1addsimp',1,S232('DSPout',dspid_left),S232('DAC',0)); %DSPout to DAC0
S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left)); %IB to DSPin
%Make connections for right ear
S232('PD1addsimp',1,S232('DSPout',dspid_right),S232('DAC',1));
S232('PD1specIB',1,S232('IB',1),S232('DSPin',dspid_right));
%Load left      
S232('pushf',ephonefilt_left,length(ephonefilt_left));
S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.ephonescale,1.0,1);
%Load right
S232('pushf',ephonefilt_right,length(ephonefilt_right));
S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.ephonescale,1.0,1);

set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4mute',1);	
S232('PA4mute',2);	

ISI = XStimParams.test_ISI;

%Get all HRTF spectra indices
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


%%%%%%%%%% load all HRTFs with HRTFfiletype == 2
dir = 0;
if XStimParams.space_flag
    tempFN = FN.space_std;
    filetype = FN.HRTFfiletype(1,2);
else
    tempFN = FN.ablequal_std;
    filetype = FN.HRTFfiletype(7,2);
end
if filetype == 2
    eval(['load -mat ' FN.space_path tempFN]);
    TF1_space = TF1; TF2_space = TF2;
    dir_space = dir;
    hrtfdirmat = dir_space;
elseif filetype == 1
    hrtfdirmat = sph2dbl(mtlrdir([FN.space_path tempFN]));
end

clear dir TF1 TF2
%%%%%%%%%%%

% find locations for 2 sources NOTE:: OFFSETS here are for loc2 vs loc1
% (loc2 = loc1 + offset)
% (not, as in Two Source, either source vs mid-array)
clear locind locations
numlocs = 0;

offset_el = XStimParams.offset_el(2);
offset_az = XStimParams.offset_az(2);
str1 = [];
for locnum = 1:size(XStimParams.locations,2)
    loc1 = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
    str1 = [str1 '  ' num2str(XStimParams.locations(1,locnum)) ' el, ' num2str(XStimParams.locations(2,locnum)) ' az, vs: ' ...
            num2str(XStimParams.locations(1,locnum) + offset_el) ' el, ' num2str(XStimParams.locations(2,locnum) + offset_az) ' az'];
    loc2 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,locnum) + offset_el) &...
        hrtfdirmat(2,:) == (XStimParams.locations(2,locnum)) + offset_az));
    if ~isempty(loc1) & ~isempty(loc2)
        numlocs = numlocs+1;
        locind(:,numlocs) = [loc1; loc2];
        locations(:,numlocs) = XStimParams.locations(:,locnum);
    end
end

if exist('locations') 
    XStimParams.locations = locations;
else
    error('no matching locations - reset and try again')
end

% write diary
if get(H.recorddata,'Value')
    if flag_reverse
        str0 = '     Mc2SrcRev  ';
        str2 = ['   SPL2:' num2str(10+XStimParams.source1_SPL) ...
            '  SPL1:' num2str(10+XStimParams.focalSPL(1)) '-' num2str(10+XStimParams.focalSPL(2))];
    else
        str0 = '     Mc2Src  ';
        str2 = ['   SPL1(' num2str(XStimParams.epoch_duration(1)) '):' num2str(10+XStimParams.source1_SPL) ...
            '  SPL2(' num2str(XStimParams.epoch_duration(2)) '):' num2str(10+XStimParams.focalSPL(1)) '-' num2str(10+XStimParams.focalSPL(2))];
    end
    tempstr = [str0 str1 str2 ...
            ' p = ' num2str(XStimParams.focalProb) '   atten: ' num2str(abs(XStimParams.curr_ABL)) ' dB'];
    update_diary
end

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

%Start clock
m110dx( C_.START);
%Send pulse: PD1 GO!
m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse


%%%%%%%%%%%%%%%%%%%% make the stimuli we'll use
set(H.status,'BackgroundColor','blue');
set(H.status,'String','Status: Building Stimuli');
set(H.remreps,'String','1');
numtrials = size(XStimParams.locations,2);
finalspikematrix = zeros(1,numtrials);

%Randomize the stimuli
randseq = randperm(numtrials);
trialnum = 1;

% make SPL sequences for each source
% source2 will jitter around in SPL, but fall to 0 SPL whenever source1 comes on
if isempty(SPLs_lo)
    XStimParams.focalProb = 1;
end

if XStimParams.repeatmode
    dur(1) = XStimParams.onEpochs + XStimParams.offEpochs;        % total epochs of each rep
    n(1) = XStimParams.onEpochs/(XStimParams.source1_ISI);        % #probes during adapter ON
    n(2) = XStimParams.offEpochs/(XStimParams.source1_ISI);       % #probes during adapter OFF
else
    dur(1) = nEpochs;        % total epochs of each rep
    n(1) = round(nEpochs/(XStimParams.source1_ISI));        % #probes during adapter ON    
end
temp_ones = ones(nEpochPts,1);

% to hold AC sequences for each location and each rep
ACseq1 = cell(XStimParams.numreps * numtrials,1);
ACseq2 = ACseq1;

% make SPLseq
randEpochs = rand(nEpochs,1);        
randEpochs2 = round(rand(nEpochs,1)*length(SPLs_hi)+.49999999999);  
randEpochs3 = round(rand(nEpochs,1)*length(SPLs_lo)+.49999999999);  

seq1 = zeros(nEpochs,1);                             % SPL for each epoch
seq2 = zeros(nEpochs,1);                             % SPL for each epoch
seqSPL1 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt
seqSPL2 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt
rampEnv = ones(nEpochPts,1);
rampEnv(1:150) = linspace(0,1,150);
rampEnv(end-149:end) = linspace(1,0,150);
if XStimParams.probeFlag(1)
    if flag_reverse
        if randEpochs(1) <= XStimParams.focalProb
            seq1(1) = SPLs_hi(randEpochs2(1));
        else
            seq1(1) = SPLs_lo(randEpochs3(1));
        end
    else
        seq1(1) = 10+XStimParams.source1_SPL;
    end
    ind0 = 1:nEpochPts;
    if flag_reverse
        seqSPL1(ind0) = seq1(1) * rampEnv;
    else
        seqSPL1(ind0) = seq1(1) * ones(size(ind0));
    end
end
for iepoch = 2:dur(1):nEpochs
    % probes
    ind1 = [];
    if XStimParams.probeFlag(1)
        ind1 = [(iepoch-1):XStimParams.source1_ISI:iepoch+XStimParams.source1_ISI*(n(1))];
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
ACseq1{1} = seq1(:)';             % not updated - same for every rep
seqSPL1 = 10.^(seqSPL1/20)-1;
ACseq2{1} = seq2(:)';
seqSPL2 = 10.^(seqSPL2/20)-1;
if length(seqSPL1) ~= nPts
    disp(['Euston, we have a problem ...']);
end    

while (exist1('H.Mc2Sourcefig') & (trialnum <= numtrials))
    %Check for pause by user
    if pause_check  return; end
    
    %%%%%%%%%%%%%%%%%%%%%%% Stim 1
    param0 = get(H.stim_type,'Value');
    if param0 ~= 9
        stim1 = get_stim(param0,XStimParams.freq(1));
    else
        stim1 = get_stim(param0,FN.stim_path,FN.stim);
    end
    
    % remove any DCoffset
    stim1 = stim1 - mom(stim1,1);
    
    % scale to AC power == 1
    stim1 = stim1/mom(stim1,2);
    
    % modulate stim1
    if ~strcmp(XStimParams.mod_type,'None')
        Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
        stim1 = stim1 .* Envelope(:)';
    end
    %Ramp stim1
    if get(H.stim_type,'Value') ~= 9
        ramp_time = 5; %ms
        [stim1] = ramp_sound(stim1,TDT.Fs,ramp_time);
    end
    
    % remove any DCoffset
    stim1 = stim1 - mom(stim1,1);
    
    %Add in the leading silent period
    stim1 =  [zeros(1,silence_len) stim1];
    
    %Add in the trailing silent period
    stim1 =  [stim1 zeros(1,silence_len2)];
    
    %Apply HRTF filtering
    if FN.HRTFfiletype(1,2) == 1
        hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locind(1,randseq(trialnum)))-1);
        hrtf_right = mtlrch([FN.space_path FN.space_std],2*locind(1,randseq(trialnum)));
    else
        hrtf_left = TF1_space(locind(1,randseq(trialnum)),:);
        hrtf_right = TF2_space(locind(1,randseq(trialnum)),:);
    end
    if XStimParams.ABLequal_flag
        temp1 = abs(fft(hrtf_left));
        temp2 = abs(fft(hrtf_right));
        ABL1 = mean(temp1(18:94) + temp2(18:94))/2;
    end
    
    stim1_L = conv(stim1,hrtf_left);
    stim1_R = conv(stim1,hrtf_right);
    
    % remove DC offset
    stim1_L =stim1_L - round(mean(stim1_L));
    stim1_R = stim1_R - round(mean(stim1_R));
    
    % these adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
    % presented through DSPs (*.eq) at 0,0 (broadband)
    stim1_L = stim1_L * TDT.hrtf_Lfactor;
    stim1_R = stim1_R * TDT.hrtf_Rfactor;
    
    %%%%% multiply by the epoch-envelopes for loc1
    lag = silence_len+100;          % +100 is an approximation of the HRTF-delay
    source1_L = zeros(size(stim1_L));               % added 3/2/09
    source1_R = zeros(size(stim1_R));
    
    source1_L(lag:lag+nPts-1) = stim1_L(lag:lag+nPts-1) .* seqSPL1;
    source1_R(lag:lag+nPts-1) = stim1_R(lag:lag+nPts-1) .* seqSPL1;
    
    % calc ACpower
    for iepoch = 1:nEpochs
        temp = source1_L(lag+(iepoch-1)*nEpochPts:silence_len+(nEpochPts-100)+(iepoch-1)*nEpochPts);
        if any(temp)    
            ACpower(1,iepoch) = std(temp);
        else
            ACpower(1,iepoch) = 0;
        end
        temp = source1_R(lag+(iepoch-1)*nEpochPts:silence_len+(nEpochPts-100)+(iepoch-1)*nEpochPts);
        if any(temp)
            ACpower(2,iepoch) = std(temp);
        else
            ACpower(2,iepoch) = 0;
        end
    end        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Stim 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    param0 = get(H.stim_type2,'Value');
    if param0 ~= 9
        stim2 = get_stim(param0,XStimParams.freq(2));
    else
        stim2 = get_stim(param0,FN.stim_path2,FN.stim2);
    end
    % remove any DCoffset
    stim2 = stim2 - mom(stim2,1);
    % scale to ACpower == 1
    stim2 = stim2/mom(stim2,2);
    
    % modulate stim2
    if ~strcmp(XStimParams.mod_type2,'None')
        Envelope = make_env(DUR, XStimParams.mod_type2, XStimParams.mod_depth(2), XStimParams.mod_freq(2), XStimParams.mod_phase(2));
        stim2 = stim2 .* Envelope(:)';
    end
    %Ramp the stimuli
    if get(H.stim_type2,'Value') ~= 9
        ramp_time = 5; %ms
        stim2 = ramp_sound(stim2,TDT.Fs,ramp_time);
    end
    % remove any DCoffset
    stim2 = stim2 - mom(stim2,1);
    
    %Add in the leading silent period
    stim2 = [zeros(1,silence_len) stim2];
    
    %Add in the trailing silent period
    stim2 = [stim2 zeros(1,silence_len2)];
    
    %Apply HRTF filtering
    if FN.HRTFfiletype(1,2) == 1
        hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locind(2,randseq(trialnum)))-1);
        hrtf_right = mtlrch([FN.space_path FN.space_std],2*locind(2,randseq(trialnum)));
    else
        hrtf_left = TF1_space(locind(2,randseq(trialnum)),:);
        hrtf_right = TF2_space(locind(2,randseq(trialnum)),:);
    end
    % re-calc for ABLequal
    if XStimParams.ABLequal_flag
        temp1 = abs(fft(hrtf_left));
        temp2 = abs(fft(hrtf_right));
        ABL2 = mean(temp1(18:94) + temp2(18:94))/2;
        factor1 = ABL1/ABL2;
    else
        factor1 = 1;
    end        
    
    stim2_L = conv(stim2,hrtf_left*factor1);
    stim2_R = conv(stim2,hrtf_right*factor1);
    
    % remove DC offset
    stim2_L = stim2_L - round(mean(stim2_L));
    stim2_R = stim2_R - round(mean(stim2_R));
    
    % adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
    % presented through DSPs (*.eq) at 0,0 (broadband)
    stim2_L = stim2_L * TDT.hrtf_Lfactor;
    stim2_R = stim2_R * TDT.hrtf_Rfactor;
    
    %%% multiply by the epoch-envelopes for loc2
    source2_L = zeros(size(stim2_L));               % added 3/2/09
    source2_R = zeros(size(stim2_R));
    source2_L(lag:lag+nPts-1) = stim2_L(lag:lag+nPts-1) .* seqSPL2;
    source2_R(lag:lag+nPts-1) = stim2_R(lag:lag+nPts-1) .* seqSPL2;
    
    % calc ACpower
    for iepoch = 1:nEpochs
        temp = source2_L(lag+(iepoch-1)*nEpochPts:silence_len+(nEpochPts-100)+(iepoch-1)*nEpochPts);
        if any(temp)    
            ACpower(3,iepoch) = std(temp);
        else
            ACpower(3,iepoch) = 0;
        end
        temp = source2_R(lag+(iepoch-1)*nEpochPts:silence_len+(nEpochPts-100)+(iepoch-1)*nEpochPts);
        if any(temp)
            ACpower(4,iepoch) = std(temp);
        else
            ACpower(4,iepoch) = 0;
        end
    end        
    
    % add two weighted stims
    trial_left = (source1_L * XStimParams.factor(1) + source2_L * XStimParams.factor(2));
    trial_right = (source1_R * XStimParams.factor(1) + source2_R * XStimParams.factor(2));
    
    % 3/7/07  NOTE: scaleFactor set from TDT.scaleFactor or 60 above
    trial_left = trial_left * scaleFactor;
    trial_right = trial_right * scaleFactor;
    
    %pad with zeros for ephonefilter
    trial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
    trial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
    
    % USED TO save stims to disk with name of loc1. NOW JUST LOAD BUFs
    if(exist1('H.Mc2Sourcefig'));
        S232('push16',trial_left,length(trial_left));
        S232('qpop16',BUF.L1);
        %fname = ['stimbuf_left_' num2str(hrtfdirmat(1,locind(1,randseq(trialnum)))) ...
        %        '_' num2str(hrtfdirmat(2,locind(1,randseq(trialnum))))];
        %evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
        %        [' ''' stimuli_dir fname ''' ']   ',0);'];
        %eval(evalstr);

        S232('push16',trial_right,length(trial_right));
        S232('qpop16',BUF.R1);
        %fname = ['stimbuf_right_' num2str(hrtfdirmat(1,locind(1,randseq(trialnum)))) ...
        %        '_' num2str(hrtfdirmat(2,locind(1,randseq(trialnum))))];
        %evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
        %        [' ''' stimuli_dir fname ''' ']   ',0);'];
        %eval(evalstr);

    end
    
    remtrials = numtrials - trialnum;
    set(H.remtrials,'String',num2str(remtrials));
    trialnum = trialnum + 1;
    pause(0);
end 										%end loop over trials

%%%%%%%%%%%%%%%%%%%%%%%%%%% finished making sounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);
pause(1)        % in case there is a click on atten reset

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
remtrials = XStimParams.numreps * numtrials;
%loop for reps
while (exist1('H.Mc2Sourcefig') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli
    randseq = randperm(numtrials);   
    trialnum = 1-firstloop;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    % loop for trials (locations) (with one extra (silent) loop)
    while (exist1('H.Mc2Sourcefig') & (trialnum <= numtrials))
        
        %Check for pause by user
        if pause_check  return; end
        
        %Wait till PD1 is finished
        %while S232('PD1status',1) usec_delay(1000); end
        usec_delay(1000);

        if(trialnum <= numtrials)
            %fname = ['stimbuf_left_' num2str(hrtfdirmat(1,locind(1,randseq(trialnum+firstloop)))) ...
            %        '_' num2str(hrtfdirmat(2,locind(1,randseq(trialnum+firstloop))))];
            %evalstr = ['S232(''disk2dama16'',BUF.L1,'  [' ''' stimuli_dir fname ''' '] ',0);'];
            %eval(evalstr);
            %fname = ['stimbuf_right_' num2str(hrtfdirmat(1,locind(1,randseq(trialnum+firstloop)))) ...
            %        '_' num2str(hrtfdirmat(2,locind(1,randseq(trialnum+firstloop))))];
            %evalstr = ['S232(''disk2dama16'',BUF.R1,' [' ''' stimuli_dir fname ''' '] ',0);'];
            %eval(evalstr);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if repnum < XStimParams.numreps     % make next stim
            % make SPLseq
            randEpochs = rand(nEpochs,1);        
            randEpochs2 = round(rand(nEpochs,1)*length(SPLs_hi)+.49999999999);  
            randEpochs3 = round(rand(nEpochs,1)*length(SPLs_lo)+.49999999999);  
            
            seq1 = zeros(nEpochs,1);                             % SPL for each epoch
            seq2 = zeros(nEpochs,1);                             % SPL for each epoch
            seqSPL1 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt
            seqSPL2 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt
            
            if XStimParams.probeFlag(1)
                if flag_reverse
                    if randEpochs(1) <= XStimParams.focalProb
                        seq1(1) = SPLs_hi(randEpochs2(1));
                    else
                        seq1(1) = SPLs_lo(randEpochs3(1));
                    end
                else
                    seq1(1) = 10+XStimParams.source1_SPL;
                end
                ind0 = 1:nEpochPts;
                if flag_reverse
                    seqSPL1(ind0) = seq1(1) * rampEnv;
                else
                    seqSPL1(ind0) = seq1(1) * ones(size(ind0));
                end
            end
            for iepoch = 2:dur(1):nEpochs
                % probes
                ind1 = [];
                if XStimParams.probeFlag(1)
                    ind1 = [(iepoch-1):XStimParams.source1_ISI:iepoch+XStimParams.source1_ISI*(n(1))];
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
            ACseq1{trialnum*(repnum-1)+2-firstloop} = seq1(:)';             % until 12/01/08 not updated - same for every rep
            seqSPL1 = 10.^(seqSPL1/20)-1;
            ACseq2{trialnum*(repnum-1)+2-firstloop} = seq2(:)';
            seqSPL2 = 10.^(seqSPL2/20)-1;
            if length(seqSPL1) ~= nPts
                disp(['Euston, we have a problem ...']);
            end 
            
            
            
            %%% multiply by the epoch-envelopes for loc1 & loc2
            source1_L = zeros(size(stim1_L));               % added 3/2/09
            source1_R = zeros(size(stim1_R));
            source2_L = zeros(size(stim2_L));
            source2_R = zeros(size(stim2_R));
            source1_L(lag:lag+nPts-1) = stim1_L(lag:lag+nPts-1) .* seqSPL1;
            source1_R(lag:lag+nPts-1) = stim1_R(lag:lag+nPts-1) .* seqSPL1;
            source2_L(lag:lag+nPts-1) = stim2_L(lag:lag+nPts-1) .* seqSPL2;
            source2_R(lag:lag+nPts-1) = stim2_R(lag:lag+nPts-1) .* seqSPL2;
            
            % add two weighted stims
            trial_left = (source1_L * XStimParams.factor(1) + source2_L * XStimParams.factor(2));
            trial_right = (source1_R * XStimParams.factor(1) + source2_R * XStimParams.factor(2));
            
            % 3/7/07  NOTE: scaleFactor set from TDT.scaleFactor or 60 above
            trial_left = trial_left * scaleFactor;
            trial_right = trial_right * scaleFactor;
            
            %pad with zeros for ephonefilter
            trial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
            trial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
            
        end             % end of making next stim
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % USED TO save stims to disk with name of loc. NOW JUST LOAD BUFs
        if (exist1('H.Mc2Sourcefig')) & repnum < XStimParams.numreps 
            S232('push16',trial_left,length(trial_left));
            S232('qpop16',BUF.L1);
            %fname = ['stimbuf_left_' num2str(hrtfdirmat(1,locind(1,randseq(trialnum+firstloop)))) ...
            %        '_' num2str(hrtfdirmat(2,locind(1,randseq(trialnum+firstloop))))];
            %evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
            %        [' ''' stimuli_dir fname ''' ']   ',0);'];
            %eval(evalstr);
            S232('push16',trial_right,length(trial_right));
            S232('qpop16',BUF.R1);
            %fname = ['stimbuf_right_' num2str(hrtfdirmat(1,locind(1,randseq(trialnum+firstloop)))) ...
            %        '_' num2str(hrtfdirmat(2,locind(1,randseq(trialnum+firstloop))))];
            %evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
            %        [' ''' stimuli_dir fname ''' ']   ',0);'];
            %eval(evalstr);
        end
       
        %Wait till PD1 is finished
        %while S232('PD1status',1) usec_delay(1000); end
        t = toc-.1;
        %if trialnum > 1
            while t <= (DUR+ISI)/1000 
                t = toc-.1;
            end
        %end        
        
        S232('PD1stop',1);
        
        %Stop the m110 and get spikes
        if(~firstloop)							% first trial just for loading
            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, DUR/2); 			% Take DUR/2 spikes max
            ind = find(spikes ~= 0); 						% Get clock events that are spikes
            spikes = spikes(ind);
            ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
            if(exist1('H.Mc2Sourcefig') & ~isempty(spikes)) 
                spikes = spikes(ind);
                spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
                EL_trial = [EL_trial;locations(1,randseq(trialnum+firstloop))* ones(size(spikes))];
                AZ_trial = [AZ_trial;locations(2,randseq(trialnum+firstloop))* ones(size(spikes))];
                repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
                Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
            end
        end
                
        %Check for pause by user
        if pause_check  return; end
        
        if(exist1('H.Mc2Sourcefig') & (remtrials>1))       %(trialnum <= numtrials))
            
            S232('seqplay',BUF.playspec1);
            S232('PD1arm',1);
            
            %Send trigger
            %Set up MII
            m100x( C_.INIT );
            m110dx( C_.INIT );
            m110dx( C_.CLOCK, mii_us_per_sample);
            m110dx( C_.MODE, M_.PST );
            
            if (remtrials>1)      %(trialnum <= numtrials)
                while toc < ISI/1000     end
                %Start clock
                m110dx( C_.START);
                %Send pulse: PD1 GO!
                m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
                tic
            end
            
        end
        
        if firstloop
            firstloop=0;
        else
            remtrials = remtrials-1;
            finalspikematrix(randseq(trialnum)) = ...
                finalspikematrix(randseq(trialnum)) + ...
                length(spikes);
        end
        
        set(H.remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end %end loop over trials
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.Mc2Sourcefig') & ~exist1('H.finalspikeratefig'))
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
    [RLFall, n1] = plot_Mc2srcA(spikes_trial, ACseq1{repnum}, RLFall, n1);            % plots the RLF
 
    
    
    %Record Data
    if(exist1('H.Mc2Sourcefig') & get(H.recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
        record_data3(XStimParams,datamatrix, [ACseq1 ACseq2], ACpower);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.remreps,'String',num2str(remreps));
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

if 0                %(exist1('H.Mc2Sourcefig') & get(H.recorddata,'Value'))
    %Plot RLF figure
    if ~exist1('H.Mc2SourceRLFfig')
        H.Mc2SourceRLFfig = figure;
    end
    figure(H.Mc2SourceRLFfig)
    hold on
    plot_Mc2src            % loads current datafile and plots the RLF
end


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
        [stim] = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
    case 9    %Stimulus from file
        if(~exist1('stim_from_file'))
            fid = fopen([param1 param2],'r');
            stim = fread(fid,inf,'float');
            fclose(fid);
        end
        if(size(stim,1) > 1)
            stim = stim';
        end
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