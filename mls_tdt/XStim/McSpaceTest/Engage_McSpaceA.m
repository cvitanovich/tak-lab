function [] = Engage_mcSpace()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

% Engage_mcSpace
%*******************************************************************************
%	The McSpace Test operation 3/15/07
% plays co-localized Dean & McAlpine-like stimuli
% also saves sequence of ACpower for each epoch
% maximum duration about 30 secs
% uses *.eq HRTF files
%*******************************************************************************
S232('trash');
S232('dropall');

% set SPL params:
scaleFactor = 60;
if ~XStimParams.HiDynamicRange
    BN = questdlg('Use HighDynamic Range?','Dynamic Range is set to low','High','Low (normal)','High');
    if ~isempty(findstr(BN,'Low'))    scaleFactor = TDT.scaleFactor;      end
end
nEpochmsecs = XStimParams.epoch_duration(1);       % duration of each epoch at one SPL
% SPLs for low prob and for hi prob
SPLs_hi = 10+(XStimParams.focalSPL(1) : XStimParams.focalSPL(2));
SPLs_lo = setdiff([10:40], SPLs_hi);
if isempty(SPLs_lo)
    XStimParams.focalProb = 1;
    set(H.focalProb,'string','1');
end

stimuli_dir = FN.temp_stim_path;
fclose all;
if ~XStimParams.buildOnly
    eval(['delete ' stimuli_dir '*.*;']);
end

% reset non-allowed stim_type to BroadBand
if isempty(find([1 2 3 6 7 8 9 10] == get(H.stim_type,'Value')))
    set(H.stim_type,'Value',8);
    disp('Unsupported stimulus type for Mcspace test; reset to BROADBAND');
end

XStimParams.curr_stimdur = str2num(get(H.mcSpace_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.mcSpace_ISI,'String'));
XStimParams.numreps = str2num(get(H.mcSpace_numreps,'String'));
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

% get stim file?
if get(H.stim_type,'Value') == 9
    if isempty(FN.stim)         %Stimulus from file
        [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select stimulus File');
    end
    D = dir([FN.stim_path FN.stim]);
    temp = (D.bytes/4) - (XStimParams.silence_lead + XStimParams.silence_trail)*(TDT.Fs/1000) - 254;
    XStimParams.curr_stimdur = round(1000 * temp / TDT.Fs);
    set(H.mcSpace_stimFN,'String',FN.stim);
end

% calc Npts_totalplay
% pts for leading silence
DUR_silence = XStimParams.silence_lead; 	%ms
silence_len = (DUR_silence * round(TDT.Fs/1000));
% pts for trailing silence
DUR_silence2 = XStimParams.silence_trail; 	%ms
silence_len2 = (DUR_silence2 * round(TDT.Fs/1000));
% duration
DUR = XStimParams.curr_stimdur;

Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen;
XStimParams.test_type = 'mcSpace FC';        str1 = 'fc';
disp('This is a FULLY CUED McSPACE test')

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
eval(['save ' FN.temp_stim_path 'XStimParams_mcSpace_' str1 ' Temp_params;'])
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
    set(H.pickerfig,'Color', [.1 .2 .8]);
    set(H.picker_error,'visible','on');
    pause;
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end
nLocs = size(XStimParams.locations,2);
clear locind
str1 = [];
for locnum = 1:nLocs
    locind(locnum) = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
    str1 = [str1 num2str(XStimParams.locations(1,locnum)) ' el, ' ...
            num2str(XStimParams.locations(2,locnum)) ' az'];
    
end

% update display
setInfo_mcSpace;

% update diary
if get(H.mcSpace_recorddata,'Value')
    str2 = ['  ' num2str(XStimParams.focalSPL(1)) '-' num2str(XStimParams.focalSPL(2)) ',  p='  num2str(XStimParams.focalProb)];
    tempstr = ['    Mcspace-test ' str1 str2 '   atten: ' num2str(abs(ABL))];
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

%%%%%  make the stimulus
set(H.mcSpace_status,'BackgroundColor','yellow');
set(H.mcSpace_status,'ForegroundColor','blue');
set(H.mcSpace_status,'String','Status: building stimulus');

switch get(H.stim_type,'Value')    
    case 1                  % tone
        stim = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur,0);
    case 2         % GammaTones
        stim = MakeGammaTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur);
    case 3                  %Use Concurrent ToneStack if specified: 10 frequencies
        tonestackfreqs = ...
            XStimParams.bandwidth(1): round((XStimParams.bandwidth(2) - XStimParams.bandwidth(1))/9):...
            XStimParams.bandwidth(2);
        stim = MakeToneStack(TDT.Fs,tonestackfreqs,XStimParams.curr_stimdur);
    case 6                %Narrowband Noise
        stim = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
    case 7                %Narrowband Noise Stack
        stim = MakeNBNoiseStack(TDT.Fs,XStimParams.curr_stimdur,XStimParams.spec_NBNoiseStackBands);
    case 8                %Broadband Noise
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
    case 9                %Stimulus from file
        fid = fopen([FN.stim_path FN.stim],'r');
        stim = fread(fid,inf,'float');
        fclose(fid);
        set(H.mcSpace_repeatmode,'value',1)
        XStimParams.repeatmode = 1;
    case 10                %Broadband Noise then flattened
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        stim = flatten(stim);
    otherwise
        set(H.stim_type,'Value',8);
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        disp('Unsupported stimulus type for space test; reset to BROADBAND');
end
nPts = length(stim);

if get(H.stim_type,'Value') ==9     nPts = nPts - silence_len - silence_len2 -254;  end

nEpochPts = round(nEpochmsecs * TDT.Fs/1000);
nEpochs = round(nPts / nEpochPts);

if get(H.stim_type,'Value')~=9
    % remove DCoffset
    trial_left = stim - mom(stim,1);
    
    % normalize ACPower==1
    trial_left = trial_left / mom(trial_left,2);
    
    % modulate
    if get(H.stim_type,'Value') ~= 9 & ~strcmp(XStimParams.mod_type,'None')
        Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
        trial_left = trial_left .* Envelope(:)';
    end

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
    trial_left =  [zeros(1,silence_len) trial_left];
    
    %Add in the trailing silent period
    trial_left =  [trial_left zeros(1,silence_len2)];
    
    %Add in the trailing silent period
    trial_left =  [trial_left zeros(1,254)];
    
    % increment testnumber
    if(exist1('H.mcSpacefig') & get(H.mcSpace_recorddata,'Value'))
        update_dataFN;
    end
    
    %%%%%%%%% make SPLseq, load buffer, while this one plays
    nPts = nEpochs * nEpochPts;
    temp_ones = ones(nEpochPts,1);
    SPLs = 10:40; nSPLs = length(SPLs);
    
    ramppts = str2num(get(H.mcSpace_ramppts,'string'));
    if ramppts  ramp = 0:(ramppts-1);   end
    
    randEpochs = rand(nEpochs,1);        
    randEpochs2 = round(rand(nEpochs,1)*length(SPLs_hi)+.49999999999);  
    randEpochs3 = round(rand(nEpochs,1)*length(SPLs_lo)+.49999999999);  
    seq = zeros(nEpochs,1); seqSPL = zeros(1,nPts);
    for iepoch = 1:nEpochs
        if randEpochs(iepoch) <= XStimParams.focalProb
            seq(iepoch) = SPLs_hi(randEpochs2(iepoch));
        else
            seq(iepoch) = SPLs_lo(randEpochs3(iepoch));
        end
        st = (1+nEpochPts*(iepoch-1)); en = (nEpochPts*iepoch);
        seqSPL(st:en) = seq(iepoch) * temp_ones;
        if ramppts
            if iepoch > 1
                fact = (seqSPL(st)-seqSPL(st-1))/(ramppts-1); 
                seqSPL(st-(ramppts/2):st+(ramppts/2)-1) = seqSPL(st-(ramppts/2)) + ramp*fact;
            end
        end
    end
    
    % dB to linear
    ACseq(1,:) = seq(:)';
    seqSPL = 10.^(seqSPL/20)-1;
    % multiply non-zero section of stim
    trial_next = trial_left;
    trial_next(silence_len+1:silence_len+nPts) = trial_left(silence_len+1:silence_len+nPts) .*seqSPL;
else
    trial_next = stim(:)';
end

% calc ACpower
for iepoch = 1:nEpochs
    temp = trial_next(silence_len+1+(iepoch-1)*nEpochPts:silence_len + iepoch*nEpochPts);
    ACpower(1,iepoch) = std(temp);
    ACvar(1,iepoch) = std(abs(hilbert(temp)));
end        
if get(H.stim_type,'Value')==9
    ACseq = floor(20*log10(ACpower/scaleFactor));               % added scaleFactor on 11/14/07
end
if XStimParams.repeatmode         % save stim to file
    filename = [FN.stim_path 'McSpace'];
    if exist(filename)==2		% file exists
        delete(filename);
    end
    [fid,message] = fopen(filename,'a+','ieee-le');
    status=fseek(fid,0,'bof');
    count = fwrite(fid, trial_next, 'float32');
    fclose(fid);
end

next_buf = BUF.L1;
S232('push16',trial_next,Npts_totalplay);
S232('qpop16',next_buf);     
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
set(H.mcSpace_status,'String','Status: playing stimuli');
set(H.mcSpace_status,'BackgroundColor','green');
set(H.mcSpace_status,'ForegroundColor','white');

while (exist1('H.mcSpacefig') & (repnum <= XStimParams.numreps))
    %Randomize the stimulus locations
    randseq = randperm(nLocs);   
    locnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    next_seq = BUF.playseq_L;

    while (exist1('H.mcSpacefig') & (locnum <= nLocs) )
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
        
        %%%%%% construct & load next stimulus envelope
        if (repnum * locnum < XStimParams.numreps * nLocs)
            % use same stim for each rep??
            if ~XStimParams.repeatmode
                % make SPLseq, load other buffer, while this one plays
                randEpochs = rand(nEpochs,1);        
                randEpochs2 = round(rand(nEpochs,1)*length(SPLs_hi)+.49999999999);  
                randEpochs3 = round(rand(nEpochs,1)*length(SPLs_lo)+.49999999999);  
                seq = zeros(nEpochs,1); seqSPL = zeros(1,nPts);
                for iepoch = 1:nEpochs
                    if randEpochs(iepoch) <= XStimParams.focalProb
                        seq(iepoch) = SPLs_hi(randEpochs2(iepoch));
                    else
                        seq(iepoch) = SPLs_lo(randEpochs3(iepoch));
                    end
                    st = (1+nEpochPts*(iepoch-1)); en = (nEpochPts*iepoch);
                    seqSPL(st:en) = seq(iepoch) * temp_ones;
                    if ramppts
                        if iepoch > 1
                            fact = (seqSPL(st)-seqSPL(st-1))/(ramppts-1); 
                            seqSPL(st-(ramppts/2):st+(ramppts/2)-1) = seqSPL(st-(ramppts/2)) + ramp*fact;
                        end
                    end
                end
                % dB to linear
                ACseq(nLocs*(repnum-1)+locnum+1,:) = seq(:)';
                seqSPL = 10.^(seqSPL/20)-1;
                % multiply non-zero section of stim
                trial_next = trial_left;
                trial_next(silence_len+1:silence_len+nPts) = trial_left(silence_len+1:silence_len+nPts) .*seqSPL;
                % calc ACpower
                for iepoch = 1:nEpochs
                    temp = trial_next(silence_len+1+(iepoch-1)*nEpochPts:silence_len + iepoch*nEpochPts);
                    ACpower(nLocs*(repnum-1)+locnum+1,iepoch) = std(temp);
                    ACvar(nLocs*(repnum-1)+locnum+1,iepoch) = std(abs(hilbert(temp)));
                end  
            else
                ACseq(nLocs*(repnum-1)+locnum+1,:) = ACseq(1,:);
                ACpower(nLocs*(repnum-1)+locnum+1,:) = ACpower(1,:);
                ACvar(nLocs*(repnum-1)+locnum+1,:) = ACvar(1,:);
            end     % end of use repeatmode
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
        if(exist1('H.mcSpacefig') & ~isempty(spikes)) 
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
        set(H.mcSpace_remtrials,'String',num2str(remtrials));
        locnum = locnum + 1;
        %pause(0);
        while t < ISI/1000
            t = toc;
            if pause_check    return; end
        end
    end %end loop over trials
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.mcSpacefig') & ~exist1('H.mcSpace_finalspikeratefig'))
        H.mcSpace_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','Space Test Spike Rate Plot',...
            'NumberTitle','off');
        H.mcSpace_spikeaxes = axes;
    end
    figure(H.mcSpace_finalspikeratefig)
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.mcSpace_spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
    colorbar
    
    %Record Data
    datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
    if (exist1('H.mcSpacefig') & get(H.mcSpace_recorddata,'Value'))
        record_data3(XStimParams,datamatrix,ACpower,ACseq,ACvar);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.mcSpace_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
    
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.mcSpace_finalspikeratefig)
set(H.mcSpace_finalspikeratefig,'Name','Final Plot for Space Test');
plotdiam1(XStimParams.locations, interimspikerate);
set(H.mcSpace_spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
colorbar
if ~isempty(locmaxspikes)
    title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
            ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))],...
        'FontSize',8);
    prompt{1} = 'Figure Handle for plotting RLF';
    prompt{2} = 'marker symbol and color';
    default{1} ='';
    default{2} = 'bo-';
    answer = inputdlg(prompt,'Plotting params for RLF',1,default);
    if isempty(answer{1})
        FigParams.handle = '';
    else
        FigParams.handle = eval(answer{1});
    end
    FigParams.marker = answer{2};
    [H0] = plotMcStim(datamatrix,ACseq,FigParams);
end




set(H.mcSpace_status,'String','Status: results');
set(H.mcSpace_status,'BackgroundColor','blue');

set(H.exitspace,'Visible','on');
set(H.resetspace,'Visible','on');

% increment test number
if(exist1('H.mcSpacefig') & get(H.mcSpace_recorddata,'Value'))
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
    case 'LP Noise'
        %B = fir1(500,param2 /(TDT.Fs/2));
        %LP_noise = FILTFILT(B, 1, rand(Npts*2,1));
        %LP_noise = LP_noise(Npts/2+1:Npts *3/2);
        %LP_noise = LP_noise - mean(LP_noise);
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
                    set(H.mcSpace_modfile,'String',[FN.mod_path FN.mod]);
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
while (exist1('H.mcSpacefig') & get(H.pausespace,'Value'))
    pause(0);
    if(~exist1('H.mcSpacefig')) return; end         
    set(H.exitspace,'Visible','on');
    set(H.resetspace,'Visible','on');
    if(exist1('H.mcSpacefig') & get(H.resetspace,'Value') == 1)
        set(H.resetspace,'Value',0);
        set(H.pausespace,'Value',0);
        Reset_space;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_mcSpace;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end