function [] = Engage_Adapt1()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

% Engage_Adapt1
%*******************************************************************************
%	The Adapt1 Test operation 3/15/07
% also saves sequence of ACpower for each epoch
% maximum duration about 30 secs
% uses *.eq HRTF files
%
%                        ...SPL2(x)...dur(2)
% stim:   SPL1...dur(1)..                   ....SPL1...dur(1)
%                                                            ....SPL2(x)...dur
%                                                            
%*******************************************************************************
% Nov 15, 2007  added capabilities to modulate epochs1 and 2 separately.

S232('trash');
S232('dropall');

% set SPL params:
randStimFlag = get(H.Adapt1_randStimFlag,'value');
scaleFactor = 60;
if ~XStimParams.HiDynamicRange
    BN = questdlg('Use HighDynamic Range?','Dynamic Range is set to low','High','Low (normal)','High');
    if ~isempty(findstr(BN,'Low'))    scaleFactor = TDT.scaleFactor;      end
end
nEpochmsecs = XStimParams.epoch_duration;       % durations for epochs1 and 2
SPLs_1 = 10 + XStimParams.SPL1;                 % SPLs for epoch1 and epoch2
SPLs_2 = 10 + XStimParams.SPL2;

% clear old stims from disk
stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);

% reset non-allowed stim_type to BroadBand
if isempty(intersect([1 2 6 8 10],get(H.stim_type,'Value')))
    set(H.stim_type,'Value',8);
    disp('Unsupported stimulus type for Adapt1 test; reset to BROADBAND');
end

XStimParams.curr_stimdur = str2num(get(H.Adapt1_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.Adapt1_ISI,'String'));
XStimParams.numreps = str2num(get(H.Adapt1_numreps,'String'));
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
S232('make',1,BUF.playseq_L);
S232('make',2,0);
S232('qpop16',BUF.playspec);

% calc Npts_totalplay
silence_lead_pts = (XStimParams.silence_lead * round(TDT.Fs/1000));
silence_mid_pts = (XStimParams.silence_mid * round(TDT.Fs/1000));
silence_trail_pts = (XStimParams.silence_trail * round(TDT.Fs/1000));
% duration
DUR = XStimParams.curr_stimdur;

Npts_totalplay = (silence_trail_pts + silence_lead_pts + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen;
XStimParams.test_type = 'Adapt1 FC';        str1 = 'fc';
disp('This is a FULLY CUED Adapt1 test')

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

dspid_left = 0; dspid_right = 1;
S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('FFF'));           % used to be '0xFFF'
S232('dropall');
%Make connections for left ear
S232('PD1addsimp',1,S232('DSPout',dspid_left),S232('DAC',0)); %DSPout to DAC0
S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left)); %IB to DSPin
%Make connections for right ear
S232('PD1addsimp',1,S232('DSPout',dspid_right),S232('DAC',1));
S232('PD1specIB',1,S232('IB',1),S232('DSPin',dspid_right));

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

ITD = XStimParams.curr_ITD;
if(abs(ITD) > 250) return; end

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4mute',1);
S232('PA4mute',2);

ISI = XStimParams.test_ISI;
% save XStimParams for this test
Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_Adapt1_' str1 ' Temp_params;'])
clear Temp_params str1

% check whether hrtf Files picked
if isempty(FN.space_eq) | ~isempty(strfind(FN.space_eq,'std'))        % FN not yet picked
    [FN.space_eq,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.eq preferred)');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',[FN.space_path FN.space_eq]);
    end
    set(H.spacefilepb,'Value',0);
end
str1 = 'FC';

FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);

%%%%%%%
if FN.HRTFfiletype(1,1) == 1
    hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space_eq]));
    for iLoc = 1:nLocs
        hrTF1(iLoc,:) = mtlrch([FN.space_path FN.space_eq],(2*iLoc-1));
        hrTF2(iLoc,:) = mtlrch([FN.space_path FN.space_eq],(2*iLoc));
    end
elseif FN.HRTFfiletype(1,1) == 2
    dir = 0;
    eval(['load -mat ' FN.space_path FN.space_eq]);
    hrtfdirmat = dir; clear dir
    hrTF1 = TF1;    hrTF2 = TF2; clear TF*
else
    disp(['HRTFfiletype incorrect'])
    return
end

%Get HRTF indices
XStimParams.locations = GUI.locations1';
if ~length(XStimParams.locations)
    B = questdlg('Push OK, pick a location. Then <CR>.','Location not picked','OK','OK','OK');
    pause;
end
nLocs = size(XStimParams.locations,2);
clear locind
str1 = [];
for locnum = 1:nLocs
    locind(locnum) = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
    str1 = [str1 num2str(XStimParams.locations(1,locnum)) 'el ,' num2str(XStimParams.locations(2,locnum)) 'az '];
end

% update display
setInfo_Adapt1;
if get(H.Adapt1_recorddata,'Value')
    tempstr = ['    Adapt1  ' str1 '   SPL1(' num2str(XStimParams.epoch_duration(1)) '):' num2str(10+XStimParams.SPL1) ...
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
%Start clock
m110dx( C_.START);
%Send pulse: PD1 GO!
tic
m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse



%%%%%  make the first stimulus
set(H.Adapt1_status,'BackgroundColor','yellow');
set(H.Adapt1_status,'ForegroundColor','blue');
set(H.Adapt1_status,'String','Status: building stimulus');
stimtype = get(H.stim_type,'Value');

switch stimtype    
    case 1                  % tone
        stim = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur,0);
    case 2                  % gammatone
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        stim = use1_ERBfilt(stim,TDT.Fs,XStimParams.curr_freq);
    case 6                %Narrowband Noise
        stim = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
    case 8                %Broadband Noise
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
    case 10
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        stim = flatten (stim);        
    otherwise
        set(H.stim_type,'Value',8);
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        disp('Unsupported stimulus type for space test; reset to BROADBAND');
end

% remove DCoffset
trial_left = stim - mom(stim,1);

% normalize ACPower==1
trial_left = trial_left / mom(trial_left,2);

% ramp
if get(H.stim_type,'Value') ~= 9
    ramp_time = 5;                          %ms
    trial_left = ramp_sound(trial_left,TDT.Fs,ramp_time);
end

% remove DCoffset
trial_left = trial_left(:)' - mom(trial_left,1);

%scale stimuli
% 3/7/07
trial_left = trial_left * scaleFactor;

%Add in the leading silent period
trial_left =  [zeros(1,silence_lead_pts) trial_left];

%Add in the trailing silent period
trial_left =  [trial_left zeros(1,silence_trail_pts)];

%Add in a trailing silent period for HRTFs
trial_left =  [trial_left zeros(1,254)];

% increment testnumber
if(exist1('H.Adapt1fig') & get(H.Adapt1_recorddata,'Value'))
    update_dataFN;
end

%%%%%%%%% make SPLseq, load buffer, while this one plays
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
% dB to linear
temp1 = 10.^(tempACseq1/20)-1;
temp2 = 10.^(tempACseq2/20)-1;
token1 = trial_left(silence_lead_pts+1:silence_lead_pts+nEpochPts(1));
token1 = token1(:);             % 1st epoch
token2 = trial_left(silence_lead_pts+1+nEpochPts(1):silence_lead_pts+nEpochPts(1)+nEpochPts(2));
token2 = token2(:);             % 2nd epoch
Factor1a = mom(token1.*ones1,2);        % should be approx the scalefactor
Factor2a = mom(token2.*ones2,2);

start = 0; iseg = 0;
for j = 1:length(SPLs_1)
    for i = 1:length(SPLs_2)
        iseg=iseg+1;
        if ~strcmp(XStimParams.mod_type,'None')     
            Envelope = make_env(nEpochPts(1)/30, XStimParams.mod_type(1,:), XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
            Envelope = Envelope(:);
            Factor1 = Factor1a/mom(token1.*Envelope,2);
            seqSPL(start+1:start+nEpochPts(1)) = temp1(j) * Envelope' *Factor1;
        else
            seqSPL(start+1:start+nEpochPts(1)) = temp1(j) *ones1;
        end        
        
        start = start + nEpochPts(1);
        if ~strcmp(XStimParams.mod_type2,'None')     
            Envelope = make_env(nEpochPts(2)/30, XStimParams.mod_type2, XStimParams.mod_depth(2), XStimParams.mod_freq(2), XStimParams.mod_phase(2));
            Envelope = Envelope(:);
            Factor2 = Factor2a/mom(token2.*Envelope,2);
            seqSPL(start+1:start+nEpochPts(2)) = temp2(i) * Envelope' *Factor2;   
        else
            seqSPL(start+1:start+nEpochPts(2)) = temp2(i) *ones2;   
        end
        start = start+nEpochPts(2);
        ACseq1(1,iseg) = tempACseq1(j);
        ACseq2(1,iseg) = tempACseq2(i);
    end
    if j<length(SPLs_1)
        seqSPL(start+1:start+silence_mid_pts) = zeros1;
        start = start+silence_mid_pts;
    end
end

nPts = length(seqSPL);

% multiply non-zero section of stim
trial_next = trial_left;
trial_next(silence_lead_pts+1:silence_lead_pts+nPts) = trial_left(silence_lead_pts+1:silence_lead_pts+nPts) .*seqSPL;

next_buf = BUF.L1;
S232('push16',trial_next,Npts_totalplay);
S232('qpop16',BUF.L1);     
%%%%%%%

% Wait for this buffer to finish playing
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
finalspikematrix = zeros(1,nLocs);
%loop for reps
repnum = 1;
set(H.Adapt1_status,'String','Status: playing stimuli');
set(H.Adapt1_status,'BackgroundColor','green');
set(H.Adapt1_status,'ForegroundColor','white');

while (exist1('H.Adapt1fig') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli
    randseq = randperm(nLocs);   
    locnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    while (exist1('H.Adapt1fig') & (locnum <= nLocs) )
        %Check for pause by user
        if pause_check    return; end
        
        % load HRTFs
        %Load left      
        S232('pushf',hrTF1(locind(randseq(locnum)),:),255);
        S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
        %Load right
        S232('pushf',hrTF2(locind(randseq(locnum)),:),255);
        S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
        
        S232('seqplay',BUF.playspec);
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
        tic
        
        %%%%%% construct & load next stimulus
        itrial = repnum * locnum; 
        if (itrial < XStimParams.numreps * nLocs)
            switch stimtype    
                case 1                  % tone
                    stim = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur,0);
                case 2                  % gammatone
                    stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
                    stim = use1_ERBfilt(stim,TDT.Fs,XStimParams.curr_freq);
                case 6                %Narrowband Noise
                    stim = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
                case 8                %Broadband Noise
                    stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
                case 10
                    stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
                    stim = flatten (stim);        
            end
            % remove DCoffset
            trial_left = stim - mom(stim,1);
            % normalize ACPower==1
            trial_left = trial_left / mom(trial_left,2);
            
            % ramp
            trial_left = ramp_sound(trial_left,TDT.Fs,ramp_time);
            %scale stimuli
            trial_left = trial_left * scaleFactor;
            %Add in the silent periods
            trial_left =  [zeros(1,silence_lead_pts) trial_left zeros(1,silence_trail_pts) zeros(1,254)];
            
            % now add SPL envelope           
            if randStimFlag
                tempACseq1 = SPLs_1(randperm(length(SPLs_1)));
                tempACseq2 = SPLs_2(randperm(length(SPLs_2)));
            else
                tempACseq1 = SPLs_1;    
                tempACseq2 = SPLs_2;    
            end
            % dB to linear
            temp1 = 10.^(tempACseq1/20)-1;
            temp2 = 10.^(tempACseq2/20)-1;
            
            start = 0; iseg = 0;
            for j = 1:length(SPLs_1)
                for i = 1:length(SPLs_2)
                    iseg=iseg+1;
                    if ~strcmp(XStimParams.mod_type,'None')     
                        Envelope = make_env(nEpochPts(1)/30, XStimParams.mod_type(1,:), XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
                        Envelope = Envelope(:);
                        Factor1 = Factor1a/mom(token1.*Envelope,2);
                        seqSPL(start+1:start+nEpochPts(1)) = temp1(j) * Envelope' *Factor1;
                    else
                        seqSPL(start+1:start+nEpochPts(1)) = temp1(j) *ones1;
                    end        
                    
                    start = start + nEpochPts(1);
                    if ~strcmp(XStimParams.mod_type2,'None')     
                        Envelope = make_env(nEpochPts(2)/30, XStimParams.mod_type2, XStimParams.mod_depth(2), XStimParams.mod_freq(2), XStimParams.mod_phase(2));
                        Envelope = Envelope(:);
                        Factor2 = Factor2a/mom(token2.*Envelope,2);
                        seqSPL(start+1:start+nEpochPts(2)) = temp2(i) * Envelope' *Factor2;   
                    else
                        seqSPL(start+1:start+nEpochPts(2)) = temp2(i) *ones2;   
                    end
                    
                    start = start + nEpochPts(2);  
                    ACseq1(itrial+1,iseg) = tempACseq1(j);
                    ACseq2(itrial+1,iseg) = tempACseq2(i);
                end
                if j<length(SPLs_1)
                    seqSPL(start+1:start+silence_mid_pts) = zeros1;
                    start = start+silence_mid_pts;
                end
            end
            
            % multiply non-zero section of stim
            trial_next = trial_left;
            trial_next(silence_lead_pts+1:silence_lead_pts+nPts) = trial_left(silence_lead_pts+1:silence_lead_pts+nPts) .*seqSPL;
            
            % check ACpower for each epoch
            start = 0; iseg = 0;
            for j = 1:length(SPLs_1)
                for i = 1:length(SPLs_2)
                    iseg=iseg+1;
                    ACpow(itrial+1,iseg) = mom(trial_next(start+1:start+nEpochPts(1)),2);
                    start = start + nEpochPts(1);
                    iseg=iseg+1;
                    ACpow(itrial+1,iseg) = mom(trial_next(start+1:start+nEpochPts(2)),2);
                    start = start + nEpochPts(2);  
                end
                if j<length(SPLs_1)
                    iseg=iseg+1;
                    ACpow(itrial+1,iseg) = mom(trial_next(start+1:start+silence_mid_pts),2);
                    start = start+silence_mid_pts;
                end
            end
            
            if next_buf == BUF.L1   
                next_buf = BUF.R1;  next_seq = BUF.playseq_R;
            else
                next_buf = BUF.L1;  next_seq = BUF.playseq_L;
            end
            S232('push16',trial_next,Npts_totalplay);
            S232('qpop16',next_buf);     
        end
        
        t = toc;
        while t < (Npts_totalplay/TDT.Fs)
            t = toc;
            if pause_check    return; end
        end
        
        %Stop the m110 and get spikes
        m110dx( C_.STOP);
        spikes = m110dx( C_.DATA, XStimParams.curr_stimdur * 2); 			    % Take (2* dur in msec) spikes max
        ind = find(spikes ~= 0); 						% Get clock events that are spikes
        spikes = spikes(ind);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
        if(exist1('H.Adapt1fig') & ~isempty(spikes)) 
            spikes = spikes(ind);
            spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
            EL_trial = [EL_trial;hrtfdirmat(1,locind(randseq(locnum)))* ones(size(spikes))];
            AZ_trial = [AZ_trial;hrtfdirmat(2,locind(randseq(locnum)))* ones(size(spikes))];
            repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
            Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
        end
        finalspikematrix(randseq(locnum)) = finalspikematrix(randseq(locnum)) + length(spikes);
        if pause_check    return; end
        
        %%%%%% reMake play specification buffer
        S232('dpush',10);
        S232('value',0);
        S232('make',0,next_seq);
        S232('make',1,next_seq);
        S232('make',2,0);
        S232('qpop16',BUF.playspec);
        %%%%%%%%%
        
        remtrials = nLocs - locnum;
        set(H.Adapt1_remtrials,'String',num2str(remtrials));
        locnum = locnum + 1;
        %pause(0);
        while t < ISI/1000
            t = toc;
            if pause_check    return; end
        end
    end %end loop over locnums
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.Adapt1fig') & ~exist1('H.Adapt1_finalspikeratefig'))
        H.Adapt1_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','Space Test Spike Rate Plot',...
            'NumberTitle','off');
        H.Adapt1_spikeaxes = axes;
    end
    figure(H.Adapt1_finalspikeratefig)
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.Adapt1_spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
    colorbar
    
    %Record Data
    datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
    if(exist1('H.Adapt1fig') & get(H.Adapt1_recorddata,'Value'))
        record_data3(XStimParams,datamatrix,ACseq1,ACseq2);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.Adapt1_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.Adapt1_finalspikeratefig)
set(H.Adapt1_finalspikeratefig,'Name','Final Plot for Space Test');
plotdiam1(XStimParams.locations, interimspikerate);
set(H.Adapt1_spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
colorbar

% plot all epochs
if ~isempty(datamatrix)
    title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
            ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))],'FontSize',8);
    H.Adapt1_epochFigure = analyze_epoch1(datamatrix,ACseq1,ACseq2);
end

set(H.Adapt1_status,'String','Status: results');
set(H.Adapt1_status,'BackgroundColor','blue');

set(H.exitAdapt1,'Visible','on');
set(H.resetAdapt1,'Visible','on');

% increment test number
if(exist1('H.Adapt1fig') & get(H.Adapt1_recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end

%%%%%%%%%%%
function [Envelope] = make_env(DUR,mod_type,param1,param2,param3)

% param1 used as follows:
%    'Sq wave' or 'tone' or 'LP noise' or 'File': mod_depth
% param2 used as follows:
%    'Sq wave' or 'tone or 'LP noise' : mod_freq
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
                    set(H.Adapt1_modfile,'String',[FN.mod_path FN.mod]);
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
while (exist1('H.Adapt1fig') & get(H.pauseAdapt1,'Value'))
    pause(0);
    if(~exist1('H.Adapt1fig')) return; end         
    set(H.exitAdapt1,'Visible','on');
    set(H.resetAdapt1,'Visible','on');
    if(exist1('H.Adapt1fig') & get(H.resetAdapt1,'Value') == 1)
        set(H.resetAdapt1,'Value',0);
        set(H.pauseAdapt1,'Value',0);
        Reset_space;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_Adapt1;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end