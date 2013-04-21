function [] = Engage_Adapt2()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

% Engage_Adapt2
%*******************************************************************************
%	The Adapt2 Test operation 3/15/07
% also saves sequence of ACpower for each epoch
% maximum duration about 30 secs
% uses *.eq HRTF files
%
%                        ...SPL2(x)...dur(2)
% stim:   SPL1...dur(1)..                   ....SPL1...dur(1)
%                                                            ....SPL2(x)...dur
%                                                            
%*******************************************************************************

waitfor(helpdlg('Masker at loc1; Probe at loc2','CHECK OFFSET!!!!'))

S232('trash');
S232('dropall');

% set SPL params:
randStimFlag = get(H.Adapt2_randStimFlag,'value');
nEpochmsecs = XStimParams.epoch_duration;       % durations for epochs1 and 2
SPLs_1 = 10 + XStimParams.SPL1;                 % SPLs for epoch1 and epoch2
SPLs_2 = 10 + XStimParams.SPL2;

% Check for scaleFactor
scaleFactor = 60;
if ~XStimParams.HiDynamicRange
    BN = questdlg('Use HighDynamic Range?','Dynamic Range is set to low','High','Low (normal)','High');
    if ~isempty(findstr(BN,'Low'))    scaleFactor = TDT.scaleFactor;      end
end

% clear old stims from disk
stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);

% reset non-allowed stim_type to BroadBand
if isempty(intersect([1 2 6 8],get(H.stim_type,'Value')))
    set(H.stim_type,'Value',8);
    disp('Unsupported stimulus type for Adapt2 test; reset to BROADBAND');
end

XStimParams.curr_stimdur = str2num(get(H.Adapt2_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.Adapt2_ISI,'String'));
XStimParams.numreps = str2num(get(H.Adapt2_numreps,'String'));
XStimParams.reset_flag = 0;

%Specify DAMA buffers
clear BUF
BUF.L1  		    = 1;        % for playing sounds
BUF.R1  		    = 2;        % for playing sounds
BUF.playseq_L		= 5;
BUF.playseq_R		= 6;
BUF.playspec		= 7;

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
S232('allot16',BUF.playspec,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L);
S232('make',1,BUF.playseq_R);
S232('make',2,0);
S232('qpop16',BUF.playspec);

% calc Npts_totalplay
silence_lead_pts = (XStimParams.silence_lead * round(TDT.Fs/1000));
silence_mid_pts = (XStimParams.silence_mid * round(TDT.Fs/1000));
silence_trail_pts = (XStimParams.silence_trail * round(TDT.Fs/1000));
% duration
DUR = XStimParams.curr_stimdur;

Npts_totalplay = (silence_trail_pts + silence_lead_pts + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen;
XStimParams.test_type = 'Adapt2 FC';        str1 = 'FC';
disp('This is a FULLY CUED Adapt2 test')

% Make Stimulus buffers
S232('dpush',Npts_totalplay);
S232('value',0);
S232('allot16',BUF.L1,Npts_totalplay);
S232('qpop16',BUF.L1);
S232('dpush',Npts_totalplay);
S232('value',0);
S232('allot16',BUF.R1,Npts_totalplay);
S232('qpop16',BUF.R1);

% set up PD1
S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,Npts_totalplay);

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
% load earphone filters
dspid_left = 0; dspid_right = 1;
S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('FFF'));
S232('dropall');
%Make connections for left ear
S232('PD1addsimp',1,S232('DSPout',dspid_left),S232('DAC',0));           %DSPout to DAC0
S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left));              %IB to DSPin
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
% save XStimParams for this test
Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_Adapt2_' str1 ' Temp_params;'])
clear Temp_params str1

%Get all HRTF spectra indices
XStimParams.locations = GUI.locations1';
if size(XStimParams.locations,2) ~= 1
    errordlg('Adapt2 requires one and only one picked location', 'SpacePicker ERROR');
    pause;
    XStimParams.locations = GUI.locations1';
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

offset_el = XStimParams.offset_el(2);
offset_az = XStimParams.offset_az(2);

locations = XStimParams.locations(:,1);
loc1 = max(find(hrtfdirmat(1,:) == locations(1,1) & hrtfdirmat(2,:) == locations(2,1)));
loc2 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,1) + offset_el) & hrtfdirmat(2,:) == (XStimParams.locations(2,1)) + offset_az));

str1 = [num2str(locations(1,1)) ' el, ' num2str(locations(2,1)) ' az, vs: ' ...
        num2str(locations(1,1) + offset_el) ' el, ' num2str(locations(2,1) + offset_az) ' az'];


if ~isempty(loc1) & ~isempty(loc2)
    XStimParams.locations = locations;
else
    errordlg('no matching locations - reset and try again', 'SpacePicker ERROR');
    pause;
end

% write diary
setInfo_Adapt2;
if get(H.Adapt2_recorddata,'Value')
   tempstr = ['    Adapt2  ' str1 '   SPL1(' num2str(XStimParams.epoch_duration(1)) '):' num2str(10+XStimParams.SPL1) ...
            '  SPL2(' num2str(XStimParams.epoch_duration(2)) '):' num2str(10+XStimParams.SPL2) ...
            '   atten: ' num2str(XStimParams.curr_ABL)];
    update_diary
end


%%%%%% start playing first quiet sound, while building stimuli
S232('seqplay',BUF.playspec);
S232('PD1arm',1);

%Set up MII
m100x( C_.INIT );
m110dx( C_.INIT );
m110dx( C_.CLOCK, mii_us_per_sample);
m110dx( C_.MODE, M_.PST );

%Start clock & send pulse to go!
m110dx( C_.START);
tic
m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse

%%%%%%%%% make SPLseq1 & 2
nEpochPts = round(nEpochmsecs * TDT.Fs/1000);
ones1 = ones(nEpochPts(1),1);
ones2 = ones(nEpochPts(2),1);
zeros1= zeros(silence_mid_pts,1);
if randStimFlag
    tempACseq1 = SPLs_1(randperm(length(SPLs_1)));
    tempACseq2 = SPLs_2(randperm(length(SPLs_2)));
else
    tempACseq2 = SPLs_2;    
    tempACseq1 = SPLs_1;    
end
start = 0; iseg = 0;
for j = 1:length(SPLs_1)
    for i = 1:length(SPLs_2)
        iseg=iseg+1;
        seqSPL1(start+1:start+nEpochPts(1)) = tempACseq1(j) *ones1;
        seqSPL2(start+1:start+nEpochPts(1)) = 0 *ones1;
        start = start + nEpochPts(1);
        seqSPL1(start+1:start+nEpochPts(2)) = 0 *ones2;
        seqSPL2(start+1:start+nEpochPts(2)) = tempACseq2(i) *ones2;
        start = start+nEpochPts(2);
        ACseq1(1,iseg) = tempACseq1(j);
        ACseq2(1,iseg) = tempACseq2(i);
    end
    if j<length(SPLs_1)
        seqSPL1(start+1:start+silence_mid_pts) = zeros1;
        seqSPL2(start+1:start+silence_mid_pts) = zeros1;
        start = start+silence_mid_pts;
    end
end
nPts = length(seqSPL1);

% dB to linear
seqSPL1 = 10.^(seqSPL1/20)-1;
seqSPL2 = 10.^(seqSPL2/20)-1;

%%%%%%%%%%%%%%%%%%%%  make two uncorrelated stimuli
set(H.Adapt2_status,'BackgroundColor','yellow');
set(H.Adapt2_status,'ForegroundColor','blue');
set(H.Adapt2_status,'String','Status: building stimulus');
set(H.Adapt2_remreps,'String','1');

switch get(H.Adapt2_stim_type(1),'Value')    
    case 1                  % tone
        stim1 = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur,0);
    case 2                  % gammatone
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        stim1 = use1_ERBfilt(stim,TDT.Fs,XStimParams.curr_freq);
    case 6                %Narrowband Noise
        stim1 = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
    case 8                %Broadband Noise
        stim1 = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
    otherwise
        set(H.stim_type,'Value',8);
        set(H.Adapt2_stim_type(1),'Value',8)
        stim1 = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        disp('Unsupported stimulus type for ADAPT2 space test; reset to BROADBAND');
end

switch get(H.Adapt2_stim_type(2),'Value')    
    case 1                  % tone
        stim2 = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur,0);
    case 2                  % gammatone
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        stim2 = use1_ERBfilt(stim,TDT.Fs,XStimParams.curr_freq);
    case 6                %Narrowband Noise
        stim2 = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
    case 8                %Broadband Noise
        stim2 = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
    otherwise
        set(H.stim_type2,'Value',8);
        set(H.Adapt2_stim_type(2),'Value',8)
        stim2 = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        disp('Unsupported stimulus type for ADAPT2 space test; reset to BROADBAND');
end
% remove DCoffset
stim1 = stim1 - mom(stim1,1);
stim2 = stim2 - mom(stim2,1);
% normalize ACPower==1
stim1 = stim1 / mom(stim1,2);
stim2 = stim2 / mom(stim2,2);
% modulate
if ~strcmp(XStimParams.mod_type,'None')
    Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
    stim1 = stim1 .* Envelope(:)';
    stim2 = stim2 .* Envelope(:)';
end

% ramp
if get(H.stim_type,'Value') ~= 9
    ramp_time = 5;                          %ms
    stim1 = ramp_sound(stim1,TDT.Fs,ramp_time);
    stim2 = ramp_sound(stim2,TDT.Fs,ramp_time);
end

% remove DCoffset
stim1 = stim1(:)' - mom(stim1,1);
stim2 = stim2(:)' - mom(stim2,1);

%Add in the leading silent period
stim1 =  [zeros(1,silence_lead_pts) stim1];
stim2 =  [zeros(1,silence_lead_pts) stim2];

%Add in the trailing silent period
stim1 =  [stim1 zeros(1,silence_trail_pts)];
stim2 =  [stim2 zeros(1,silence_trail_pts)];

%%%%%%%%%%%%%%%%%%%  stim1 from loc1
%Apply HRTF filtering
if FN.HRTFfiletype(1,2) == 1
    hrtf_left1 = mtlrch([FN.space_path FN.space_std],(2*loc1)-1);
    hrtf_right1 = mtlrch([FN.space_path FN.space_std],2*loc1);
else
    hrtf_left1 = TF1_space(loc1,:);
    hrtf_right1 = TF2_space(loc1,:);
end
if XStimParams.ABLequal_flag
    temp1 = abs(fft(hrtf_left1));
    temp2 = abs(fft(hrtf_right1));
    ABL1 = mean(temp1(18:94) + temp2(18:94))/2;
end
source1_L = conv(stim1,hrtf_left1);
source1_R = conv(stim1,hrtf_right1);

% remove DC offset
source1_L = source1_L - round(mean(source1_L));
source1_R = source1_R - round(mean(source1_R));

% these adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
% presented through DSPs (*.eq) at 0,0 (broadband)
source1_L = source1_L * TDT.hrtf_Lfactor;
source1_R = source1_R * TDT.hrtf_Rfactor;

%%%%% multiply by the epoch-envelopes for loc1 (probe)
lag = silence_lead_pts+100;          % +100 is an approximation of the HRTF-delay
source1_L(lag:lag+nPts-1) = source1_L(lag:lag+nPts-1) .* seqSPL1;
source1_R(lag:lag+nPts-1) = source1_R(lag:lag+nPts-1) .* seqSPL1;

%%%%%%%%%%%%%%%%%%%  stim2 from loc2
%Apply HRTF filtering
if FN.HRTFfiletype(1,2) == 1
    hrtf_left2 = mtlrch([FN.space_path FN.space_std],(2*loc2)-1);
    hrtf_right2 = mtlrch([FN.space_path FN.space_std],2*loc2);
else
    hrtf_left2 = TF1_space(loc2,:);
    hrtf_right2 = TF2_space(loc2,:);
end
% re-calc for ABLequal
if XStimParams.ABLequal_flag
    temp1 = abs(fft(hrtf_left2));
    temp2 = abs(fft(hrtf_right2));
    ABL2 = mean(temp1(18:94) + temp2(18:94))/2;
    factor1 = ABL1/ABL2;
else
    factor1 = 1;
end        
source2_L = conv(stim2,hrtf_left2*factor1);
source2_R = conv(stim2,hrtf_right2*factor1);

% remove DC offset
source2_L = source2_L - round(mean(source2_L));
source2_R = source2_R - round(mean(source2_R));

% these adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
% presented through DSPs (*.eq) at 0,0 (broadband)
source2_L = source2_L * TDT.hrtf_Lfactor;
source2_R = source2_R * TDT.hrtf_Rfactor;

%%%%% multiply by the epoch-envelopes for loc2 (adapter)
lag = silence_lead_pts+100;          % +100 is an approximation of the HRTF-delay
source2_L(lag:lag+nPts-1) = source2_L(lag:lag+nPts-1) .* seqSPL2;
source2_R(lag:lag+nPts-1) = source2_R(lag:lag+nPts-1) .* seqSPL2;

%%%%%%%%%%%%%%%%%%%    % add two weighted stims
trial_left = (source1_L * XStimParams.factor(1) + source2_L * XStimParams.factor(2));
trial_right = (source1_R * XStimParams.factor(1) + source2_R * XStimParams.factor(2));

% 3/7/07  NOTE: scaleFactor set from TDT.scaleFactor or 60 above
trial_left = trial_left * scaleFactor;
trial_right = trial_right * scaleFactor;

%pad with zeros for ephonefilter
trial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
trial_right = [trial_right zeros(1,TDT.ephonefiltlen)];

% load stims to DAMAs
S232('push16',trial_left,length(trial_left));
S232('qpop16',BUF.L1);
S232('push16',trial_right,length(trial_right));
S232('qpop16',BUF.R1);
%%%%%%%%%%%%%% finished making stmuli

% update testnumber
if(exist1('H.Mc2Sourcefig') & get(H.recorddata,'Value'))
    update_dataFN;
end

% Wait for this buffer to finish playing, then stop PD1 and read spikes
t = toc;
while t < (Npts_totalplay/TDT.Fs)+ISI/1000
    t = toc;
end
m110dx( C_.STOP);
spikes = m110dx( C_.DATA, 100); 			% Take 100 spikes max

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);
pause(1)        % in case reset atten gives click

datamatrix = [];
finalspikematrix = zeros(1,1);
%loop for reps
repnum = 1;
set(H.Adapt2_status,'String','Status: playing stimuli');
set(H.Adapt2_status,'BackgroundColor','green');
set(H.Adapt2_status,'ForegroundColor','white');

while (exist1('H.Adapt2fig') & (repnum <= XStimParams.numreps))
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    %Check for pause by user
    if pause_check    return; end
    
    S232('seqplay',BUF.playspec);
    S232('PD1arm',1);
    
    %Send trigger
    %Set up MII
    m100x( C_.INIT );
    m110dx( C_.INIT );
    m110dx( C_.CLOCK, mii_us_per_sample);
    m110dx( C_.MODE, M_.PST );
    
    %Start clock & send pulse to go!
    m110dx( C_.START);
    m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
    tic
    
    %%%%%% NOW PLAYING STIM
    if repnum < XStimParams.numreps
        
        %%%%%%%%%%%%%%%%%% construct next stimulus 
                %%%%%%%%% make SPLseq1 & 2
        nEpochPts = round(nEpochmsecs * TDT.Fs/1000);
        ones1 = ones(nEpochPts(1),1);
        ones2 = ones(nEpochPts(2),1);
        zeros1= zeros(silence_mid_pts,1);
        if randStimFlag
            tempACseq1 = SPLs_1(randperm(length(SPLs_1)));
            tempACseq2 = SPLs_2(randperm(length(SPLs_2)));
        else
            tempACseq2 = SPLs_2;    
            tempACseq1 = SPLs_1;    
        end
        start = 0; iseg = 0;
        for j = 1:length(SPLs_1)
            for i = 1:length(SPLs_2)
                iseg=iseg+1;
                seqSPL1(start+1:start+nEpochPts(1)) = tempACseq1(j) *ones1;
                seqSPL2(start+1:start+nEpochPts(1)) = 0 *ones1;
                start = start + nEpochPts(1);
                seqSPL1(start+1:start+nEpochPts(2)) = 0 *ones2;
                seqSPL2(start+1:start+nEpochPts(2)) = tempACseq2(i) *ones2;
                start = start+nEpochPts(2);
                ACseq1(repnum+1,iseg) = tempACseq1(j);
                ACseq2(repnum+1,iseg) = tempACseq2(i);
            end
            if j<length(SPLs_1)
                seqSPL1(start+1:start+silence_mid_pts) = zeros1;
                seqSPL2(start+1:start+silence_mid_pts) = zeros1;
                start = start+silence_mid_pts;
            end
        end

% dB to linear
seqSPL1 = 10.^(seqSPL1/20)-1;
seqSPL2 = 10.^(seqSPL2/20)-1;

% make two uncorrelated stimuli        
switch get(H.Adapt2_stim_type(1),'Value')    
    case 1                  % tone
        stim1 = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur,0);
    case 2                  % gammatone
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        stim1 = use1_ERBfilt(stim,TDT.Fs,XStimParams.curr_freq);
    case 6                %Narrowband Noise
        stim1 = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
    case 8                %Broadband Noise
        stim1 = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
end

switch get(H.Adapt2_stim_type(2),'Value')    
    case 1                  % tone
        stim2 = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur,0);
    case 2                  % gammatone
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        stim2 = use1_ERBfilt(stim,TDT.Fs,XStimParams.curr_freq);
    case 6                %Narrowband Noise
        stim2 = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
    case 8                %Broadband Noise
        stim2 = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
end

% remove DCoffset
stim1 = stim1 - mom(stim1,1);
stim2 = stim2 - mom(stim2,1);
% normalize ACPower==1
stim1 = stim1 / mom(stim1,2);
stim2 = stim2 / mom(stim2,2);
% modulate
if ~strcmp(XStimParams.mod_type,'None')
    Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
    stim1 = stim1 .* Envelope(:)';
    stim2 = stim2 .* Envelope(:)';
end

% ramp
if get(H.stim_type,'Value') ~= 9
    ramp_time = 5;                          %ms
    stim1 = ramp_sound(stim1,TDT.Fs,ramp_time);
    stim2 = ramp_sound(stim2,TDT.Fs,ramp_time);
end

% remove DCoffset
stim1 = stim1(:)' - mom(stim1,1);
stim2 = stim2(:)' - mom(stim2,1);

%Add in the leading silent period
stim1 =  [zeros(1,silence_lead_pts) stim1];
stim2 =  [zeros(1,silence_lead_pts) stim2];

%Add in the trailing silent period
stim1 =  [stim1 zeros(1,silence_trail_pts)];
stim2 =  [stim2 zeros(1,silence_trail_pts)];

%%%%%%%%%%%%%%%%%%%  stim1 from loc1
%Apply HRTF filtering
source1_L = conv(stim1,hrtf_left1);
source1_R = conv(stim1,hrtf_right1);

% remove DC offset
source1_L = source1_L - round(mean(source1_L));
source1_R = source1_R - round(mean(source1_R));

% these adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
% presented through DSPs (*.eq) at 0,0 (broadband)
source1_L = source1_L * TDT.hrtf_Lfactor;
source1_R = source1_R * TDT.hrtf_Rfactor;

%%%%% multiply by the epoch-envelopes for loc1 (adapter)
lag = silence_lead_pts+100;          % +100 is an approximation of the HRTF-delay
source1_L(lag:lag+nPts-1) = source1_L(lag:lag+nPts-1) .* seqSPL1;
source1_R(lag:lag+nPts-1) = source1_R(lag:lag+nPts-1) .* seqSPL1;

%%%%%%%%%%%%%%%%%%%  stim2 from loc2
%Apply HRTF filtering
source2_L = conv(stim2,hrtf_left2*factor1);
source2_R = conv(stim2,hrtf_right2*factor1);

% remove DC offset
source2_L = source2_L - round(mean(source2_L));
source2_R = source2_R - round(mean(source2_R));

% these adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
% presented through DSPs (*.eq) at 0,0 (broadband)
source2_L = source2_L * TDT.hrtf_Lfactor;
source2_R = source2_R * TDT.hrtf_Rfactor;

%%%%% multiply by the epoch-envelopes for loc2 (probe)
lag = silence_lead_pts+100;          % +100 is an approximation of the HRTF-delay
source2_L(lag:lag+nPts-1) = source2_L(lag:lag+nPts-1) .* seqSPL2;
source2_R(lag:lag+nPts-1) = source2_R(lag:lag+nPts-1) .* seqSPL2;

%%%%%%%%%%%%%%%%%%%    % add two weighted stims
trial_left = (source1_L * XStimParams.factor(1) + source2_L * XStimParams.factor(2));
trial_right = (source1_R * XStimParams.factor(1) + source2_R * XStimParams.factor(2));

% 3/7/07  NOTE: scaleFactor set from TDT.scaleFactor or 60 above
trial_left = trial_left * scaleFactor;
trial_right = trial_right * scaleFactor;

%pad with zeros for ephonefilter
trial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
trial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
%%%%%%%%%%%%%% finished making stmuli

end
    
    %Wait till PD1 is finished
    t = toc-.1;
    while t <= DUR/1000 t = toc-.1;    end
    
    S232('PD1stop',1);
    m110dx( C_.STOP);
    %Stop the m110 and get spikes
    
    % load next stim
    if repnum < XStimParams.numreps
        S232('push16',trial_left,length(trial_left));
        S232('qpop16',BUF.L1);
        S232('push16',trial_right,length(trial_right));
        S232('qpop16',BUF.R1);
    end
    
    % read spikes
    spikes = m110dx( C_.DATA, DUR/2); 			% Take DUR/2 spikes max
    ind = find(spikes ~= 0); 						% Get clock events that are spikes
    spikes = spikes(ind);
    ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
    if(exist1('H.Adapt2fig') & ~isempty(spikes)) 
        spikes = spikes(ind);
        spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
        EL_trial = [EL_trial;locations(1,1)* ones(size(spikes))];
        AZ_trial = [AZ_trial;locations(2,1)* ones(size(spikes))];
        repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
        Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
    end
    
    %Check for pause by user
    if pause_check  return; end
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.Adapt2fig') & ~exist1('H.Adapt2_finalspikeratefig'))
        H.Adapt2_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','Space Test Spike Rate Plot',...
            'NumberTitle','off');
        H.Adapt2_spikeaxes = axes;
    end
    figure(H.Adapt2_finalspikeratefig)
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.Adapt2_spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
    colorbar
    
    %Record Data
    datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
    if(exist1('H.Adapt2fig') & get(H.Adapt2_recorddata,'Value'))
        record_data3(XStimParams,datamatrix,ACseq1,ACseq2);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.Adapt2_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.Adapt2_finalspikeratefig)
set(H.Adapt2_finalspikeratefig,'Name','Final Plot for Space Test');
plotdiam1(XStimParams.locations, interimspikerate);
set(H.Adapt2_spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
colorbar

% plot all epochs
if ~isempty(datamatrix)
    title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
            ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))],'FontSize',8);
    H.Adapt2_epochFigure = analyze_epoch1(datamatrix,ACseq1,ACseq2);
end

set(H.Adapt2_status,'String','Status: results');
set(H.Adapt2_status,'BackgroundColor','blue');

set(H.exitAdapt2,'Visible','on');
set(H.resetAdapt2,'Visible','on');

% increment test number
if(exist1('H.Adapt2fig') & get(H.Adapt2_recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end


%%%%%%%%%%%
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
    case 'Sq wave'
        T = 0:Inc:(DUR/1000 - Inc);
        Tone = (param1 / 2)* square(2 * pi * param2 .* T + (.75 * 2 * pi + param3));
        Envelope = Tone + (1-param1/2);
    case 'Tone'
        T = 0:Inc:(DUR/1000 - Inc);
        Tone = (param1 / 2)* sin(2 * pi * param2 .* T + (.75 * 2 * pi + param3));
        Envelope = Tone + (1-param1/2);
    case 'LP Noise'	                            % not functional yet
        B = fir1(500,param2 /(TDT.Fs/2));
        LP_noise = FILTFILT(B, 1, rand(Npts*2,1));
        LP_noise = LP_noise(Npts/2+1:Npts *3/2);
        LP_noise = LP_noise - mean(LP_noise);
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
                    set(H.Adapt2_modfile,'String',[FN.mod_path FN.mod]);
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
%%%%%%%%%
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
while (exist1('H.Adapt2fig') & get(H.pauseAdapt2,'Value'))
    pause(0);
    if(~exist1('H.Adapt2fig')) return; end         
    set(H.exitAdapt2,'Visible','on');
    set(H.resetAdapt2,'Visible','on');
    if(exist1('H.Adapt2fig') & get(H.resetAdapt2,'Value') == 1)
        set(H.resetAdapt2,'Value',0);
        set(H.pauseAdapt2,'Value',0);
        Reset_space;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_Adapt2;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end