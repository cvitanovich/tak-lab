function [] = Engage_Space2()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

% Engage_Space2
%*******************************************************************************
%	The Space Test operation 3/15/07
% eliminated double buffering
% does not write to files
% maximum duration about 30 secs
%*******************************************************************************

stimuli_dir = FN.temp_stim_path;
fclose all;
if ~XStimParams.buildOnly
    eval(['delete ' stimuli_dir '*.*;']);
end

% reset non-allowed stim_type to BroadBand
if isempty(find([1 2 3 6 7 8 9] == get(H.stim_type,'Value')))
    set(H.stim_type,'Value',8);
    disp('Unsupported stimulus type for space test; reset to BROADBAND');
end

set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

XStimParams.curr_ITD = str2num(get(H.space2_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.space2_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.space2_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space2_ISI,'String'));
XStimParams.numreps = str2num(get(H.space2_numreps,'String'));
XStimParams.reset_flag = 0;

% load filter coeffs & other params for ABLalone test
if XStimParams.ABLalone_flag
    Fs = 30000;
    cF = round(1000*exp(([40:168]/48)*log(2)))'; 
    n_cF = length(cF);
    fcoefs = Make_ERBFiltA(Fs,cF);
    maxFactor = .00003764*cF(n_cF)+.6236;
    Factor = maxFactor ./ (.00003764*cF+.6236);
    Factormat = repmat(Factor,1,255);
    fftpts = 2048;
    freq = 0:15000/((fftpts/2)-1):15000;
    freq_ind = nearest_index(freq,cF);
    endpt = max1(freq_ind);
    startpt = min1(freq_ind);
    Xpart1 = startpt:endpt;
    Xpart2 = (fftpts + 2 - endpt):(fftpts + 2 - startpt);
end

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

% get stim file?
if get(H.stim_type,'Value') == 9
    if isempty(FN.stim)         %Stimulus from file
        [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select stimulus File');
    end
    D = dir([FN.stim_path FN.stim]);
    XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
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
if XStimParams.ildalone_flag 				%ILDAlone Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen + TDT.itdfiltlen;
    XStimParams.test_type = 'Space2 ILA';       str1 = 'ila';
    disp('This is an ILD-ALONE SPACE test')
elseif XStimParams.itdalone_flag 			%ITDAlone Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + (2* TDT.hrtffiltlen);
    XStimParams.test_type = 'Space2 ITA';       str1 = 'ita';
    disp('This is an ITD-ALONE SPACE test')
elseif XStimParams.space_flag 			%fully-cued Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen;
    XStimParams.test_type = 'Space2 FC';        str1 = 'fc';
    disp('This is a FULLY CUED SPACE test')
elseif XStimParams.ABLalone_flag			%ABL-alone Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen;
    XStimParams.test_type = 'Space2 ABA';       str1 = 'aba';
    disp('This is an ABL-alone SPACE test')
end   

% Make Stimulus buffers
S232('allot16',BUF.L1,Npts_totalplay);
S232('allot16',BUF.R1,Npts_totalplay);

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
if get(H.space2_recorddata,'Value')
    tempstr = [str1 ' space2-test: ' num2str(length(XStimParams.locations)) ' locations;     atten: ' num2str(abs(ABL))];
    update_diary
end
Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_Space2_' str1 ' Temp_params;'])
clear Temp_params str1

% check whether hrtf Files picked
if XStimParams.space_flag | XStimParams.ABLalone_flag
    if isempty(FN.space_std) | ~isempty(strfind(FN.space_std,'eq'))        % FN not yet picked
        [FN.space_std,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.std preferred)');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space_std]);
        end
        set(H.spacefilepb,'Value',0);
    end
    t1 = 1;   tempPATH = FN.space_path;   tempFN = FN.space_std;    str1 = 'FC';
elseif XStimParams.ildalone_flag
    if isempty(FN.ildalone_std)        % FN not yet picked
        [FN.ildalone_std,FN.ILA_path] = uigetfile('*.*','Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.space2_ildalonefile,'String',[FN.ILA_path FN.ildalone_std]);
        end
        set(H.space2_ildalonepb,'Value',0);
    end
    t1 = 2;   tempPATH = FN.ILA_path;   tempFN = FN.ildalone_std;    str1 = 'ILA';
elseif XStimParams.itdalone_flag        % need both, ila and ita
    if isempty(FN.ildalone_std)        % ila FN not yet picked
        [FN.ildalone_std,FN.ILA_path] = uigetfile('*.*','Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.space2_ildalonefile,'String',[FN.ILA_path FN.ildalone_std]);
        end
    end
    if isempty(FN.itdalone_std)        % ita FN not yet picked
        [FN.itdalone_std,FN.ITA_path] = uigetfile('*.*','Select ITD-alone HRTF File');
        if(FN.ITA_path ~= 0)
            set(H.space2_itdalonefile,'String',[FN.ITA_path FN.itdalone_std]);
        end
        set(H.space2_itdalonepb,'Value',0);
    end
    t1 = 3;   tempPATH = FN.ITA_path;   tempFN = FN.itdalone_std;    str1 = 'ITA';
end   
FN.HRTFfiletype(t1,2) = testHRTFfiletype(tempPATH, tempFN);

%%%%%%%
if FN.HRTFfiletype(t1,2) == 1
    hrtfdirmat = sph2dbl(mtlrdir([tempPATH tempFN]));
    for iLoc = 1:nLocs
        hrTF1(iLoc,:) = mtlrch([tempPATH tempFN],(2*iLoc-1));
        hrTF2(iLoc,:) = mtlrch([tempPATH tempFN],(2*iLoc));
    end
elseif FN.HRTFfiletype(t1,2) == 2
    dir = 0;
    eval(['load -mat ' tempPATH tempFN]);
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
for locnum = 1:nLocs
    locind(locnum) = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
end





% modify HRTFs for ABLalone test to ABL at ref loc
if XStimParams.ABLalone_flag
    ind = max(find(hrtfdirmat(1,:) == XStimParams.el & hrtfdirmat(2,:) == XStimParams.az));
    tempL = ERBFilterBankB(hrTF1(ind,:)', fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
    tempR = ERBFilterBankB(hrTF2(ind,:)', fcoefs) .* Factormat;
    [temp ablref] = calclevel_time(tempL,tempR, cF);
    FT_L = fft(hrTF1(ind,:)',fftpts);
    FT_R = fft(hrTF2(ind,:)',fftpts);
    for iLoc = 1:nLocs
        tempL = ERBFilterBankB(hrTF1(iLoc,:)', fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
        tempR = ERBFilterBankB(hrTF2(iLoc,:)', fcoefs) .* Factormat;
        [ildx ablx] = calclevel_time(tempL,tempR, cF);
        ABLfactor = (10 .^(ablx/20)) ./ (10 .^(ablref/20));
        ABLfactor_long = (interp1(freq(freq_ind),ABLfactor,freq(startpt:endpt)))';
        % apply ABLfactor to reference location
        FT_L(Xpart1)= FT_L(Xpart1) .* (ABLfactor_long);       	% positive freqs
        FT_L(Xpart2)= FT_L(Xpart2) .* flipud(ABLfactor_long);   % negative freqs
        temp = real(ifft(FT_L));      hrTF1(iLoc,:) = temp(1:255);
        FT_R(Xpart1)= FT_R(Xpart1) .* (ABLfactor_long);       	% positive freqs
        FT_R(Xpart2)= FT_R(Xpart2) .* flipud(ABLfactor_long);   % negative freqs
        temp = real(ifft(FT_R));      hrTF2(iLoc,:) = temp(1:255);
    end
    clear tempL tempR temp
end

% update display
setInfo_space2;

%%%%%  make the stimulus
set(H.space2_status,'BackgroundColor','yellow');
set(H.space2_status,'ForegroundColor','blue');
set(H.space2_status,'String','Status: building stimulus');


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
    otherwise
        set(H.stim_type,'Value',8);
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        disp('Unsupported stimulus type for space test; reset to BROADBAND');
end

% remove DCoffset
trial_left = stim - mom(stim,1);

% normalize ACPower==1
trial_left = trial_left / mom(trial_left,2);

% modulate
if get(H.stim_type,'Value') ~= 9 & ~strcmp(XStimParams.mod_type,'None')
    Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
    trial_left = trial_left .* Envelope(:)';
end

% ramp
if get(H.stim_type,'Value') ~= 9
    ramp_time = 5;                          %ms
    trial_left = ramp_sound(trial_left,TDT.Fs,ramp_time);
end

% remove DCoffset
trial_left = trial_left - mom(trial_left,1);
% make right sound
trial_right =trial_left;

% Apply ITD filtering if conducting ILDAlone Space Test
if (XStimParams.ildalone_flag)
    itdleft = 0; itdright = 0;
    ITD = round(str2num(get(H.space2_ITD,'String')));
    if(ITD < 0)
        itdleft = 0;     itdright = abs(ITD);
    else
        itdleft = abs(ITD);     itdright = 0;
    end
    eval(['load ' FN.ITD_path 'itdfilt' num2str(itdleft)]);
    eval(['filt_left = itd_filt' num2str(itdleft) ';']);
    filt_left(1) = 0;   filt_left(end) = 0;
    eval(['load ' FN.ITD_path 'itdfilt' num2str(itdright)]);
    eval(['filt_right = itd_filt' num2str(itdright) ';']);
    filt_right(1) = 0;  filt_right(end) = 0;
    trial_left = conv(trial_left,filt_left);
    trial_right = conv(trial_right,filt_right);
end

% Apply ILD filtering if conducting ITDalone Test
if XStimParams.itdalone_flag
    if FN.HRTFfiletype(2,2) == 1
        dir = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone_std]));
        ind = max(find(dir(1,:) == XStimParams.el & dir(2,:) == XStimParams.az));
        filt_left = mtlrch([FN.ILA_path FN.ildalone_std],(2*ind-1));
        filt_right = mtlrch([FN.ILA_path FN.ildalone_std],(2*ind));
    else
        dir = 0;
        eval(['load -mat ' FN.ILA_path FN.ildalone_std]);
        ind = max(find(dir(1,:) == XStimParams.el & dir(2,:) == XStimParams.az));
        filt_left = TF1(ind,:);    filt_right = TF2(ind,:); clear dir TF*
    end
    filt_left(1) = 0;   filt_left(end) = 0;
    filt_right(1) = 0;   filt_right(end) = 0;
    trial_left = conv(trial_left,filt_left);
    trial_right = conv(trial_right,filt_right);
end

% remove DC offset
trial_left = trial_left - round(mean(trial_left));
trial_right = trial_right - round(mean(trial_right));

%scale stimuli
% 3/7/07
trial_left = trial_left * TDT.scaleFactor;
trial_right = trial_right * TDT.scaleFactor;

%Add in the leading silent period
trial_left =  [zeros(1,silence_len) trial_left];
trial_right = [zeros(1,silence_len) trial_right];

%Add in the trailing silent period
trial_left =  [trial_left zeros(1,silence_len2)];
trial_right = [trial_right zeros(1,silence_len2)];

%Add in the trailing silent period
trial_left =  [trial_left zeros(1,254)];
trial_right = [trial_right zeros(1,254)];

% load to buffers
S232('push16',trial_left,Npts_totalplay);
S232('qpop16',BUF.L1);
S232('push16',trial_right,Npts_totalplay);
S232('qpop16',BUF.R1);

%Begin playing sounds   
set(H.space2_status,'String','Status: playing silent PRE_stimulus');
set(H.space2_status,'BackgroundColor','green');
set(H.space2_status,'ForegroundColor','white');
%set(H.space2_remreps,'String',num2str(XStimParams.numreps));

% increment testnumber
if(exist1('H.space2fig') & get(H.space2_recorddata,'Value'))
    update_dataFN;
end

%%%%% play out one sound to get rid of duration-change artifact
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

% Wait for these buffers to finish playing
t = toc;
while t < (Npts_totalplay/TDT.Fs)+ISI/1000
    t = toc;
end

m110dx( C_.STOP);
spikes = m110dx( C_.DATA, 100); 			% Take 100 spikes max
%%%%%%%%%%%%%%%%%%%%%%%% main data collection
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

datamatrix = [];
finalspikematrix = zeros(1,nLocs);
%loop for reps
repnum = 1;
set(H.space2_status,'String','Status: playing stimuli');
set(H.space2_status,'BackgroundColor','green');
set(H.space2_status,'ForegroundColor','white');

while (exist1('H.space2fig') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli
    randseq = randperm(nLocs);   
    trialnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    while (exist1('H.space2fig') & (trialnum <= nLocs) )
        %Check for pause by user
        if pause_check    return; end
        
        % load HRTFs
        %Load left      
        S232('pushf',hrTF1(locind(randseq(trialnum)),:),255);
        S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
        %Load right
        S232('pushf',hrTF2(locind(randseq(trialnum)),:),255);
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
        
        t = toc;
        while t < (Npts_totalplay/TDT.Fs)           %+ISI/1000
            t = toc;
            if pause_check    return; end
        end
        t = toc;
        %Stop the m110 and get spikes
        m110dx( C_.STOP);
        spikes = m110dx( C_.DATA, XStimParams.curr_stimdur * 2); 			    % Take (2* dur in msec) spikes max
        ind = find(spikes ~= 0); 						% Get clock events that are spikes
        spikes = spikes(ind);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
        if(exist1('H.space2fig') & ~isempty(spikes)) 
            spikes = spikes(ind);
            spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
            EL_trial = [EL_trial;hrtfdirmat(1,locind(randseq(trialnum)))* ones(size(spikes))];
            AZ_trial = [AZ_trial;hrtfdirmat(2,locind(randseq(trialnum)))* ones(size(spikes))];
            repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
            Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
        end
        finalspikematrix(randseq(trialnum)) = finalspikematrix(randseq(trialnum)) + length(spikes);
        if pause_check    return; end
        
        remtrials = nLocs - trialnum;
        set(H.space2_remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        %pause(0);
        while t < ISI/1000
            t = toc;
            if pause_check    return; end
        end
    end %end loop over trials
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.space2fig') & ~exist1('H.space2_finalspikeratefig'))
        H.space2_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','Space Test Spike Rate Plot',...
            'NumberTitle','off');
        H.space2_spikeaxes = axes;
    end
    figure(H.space2_finalspikeratefig)
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.space2_spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
    colorbar
    
    %Record Data
    if(exist1('H.space2fig') & get(H.space2_recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
        record_data3(XStimParams,datamatrix);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.space2_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
    
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.space2_finalspikeratefig)
set(H.space2_finalspikeratefig,'Name','Final Plot for Space Test');
plotdiam1(XStimParams.locations, interimspikerate);
set(H.space2_spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
        ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))],...
    'FontSize',8);
colorbar

% plot weights
if get(H.space2_plotwts,'Value')
    set(H.space2_status,'String','Status: calculating weights');
    set(H.space2_status,'BackgroundColor','red');
    plot_MLSweights(finalspikematrix,locind);
    set(H.space2_plotwts,'value',0);
end

set(H.space2_status,'String','Status: results');
set(H.space2_status,'BackgroundColor','blue');

set(H.exitspace,'Visible','on');
set(H.resetspace,'Visible','on');

% increment test number
if(exist1('H.space2fig') & get(H.space2_recorddata,'Value'))
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
                    set(H.space2_modfile,'String',[FN.mod_path FN.mod]);
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
while (exist1('H.space2fig') & get(H.pausespace,'Value'))
    pause(0);
    if(~exist1('H.space2fig')) return; end         
    set(H.exitspace,'Visible','on');
    set(H.resetspace,'Visible','on');
    if(exist1('H.space2fig') & get(H.resetspace,'Value') == 1)
        set(H.resetspace,'Value',0);
        set(H.pausespace,'Value',0);
        Reset_space;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_space2;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end