function [] = Engage_space4()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

% Engage_space4
%*******************************************************************************
%	The space4 Test operation 3/15/07
% also saves stimulus for each epoch
% maximum duration about 30 secs
% uses *.eq HRTF files
% ONLY ALLOWS ONE LOCATION and only FULLY CUED
% at 'focalProb' the same sound is repeated (instead of created anew)
% param3 provides the rep# for each sound (rep#1 is the repeat)
% mu sets the mean level in dB for each epoch (in conjunction with atten)
% sigma sets the peak-to-trough modulation in dB of the envelope
% before exponentiation (to dB), the envelope has a nearly uniform
% amplitude distribution (thanks to Monte Escabi)
% trying to give Nagel&Doupe-like stimuli
%*******************************************************************************

S232('trash');
S232('dropall');

% set scaleFactor
scaleFactor = 60;
if ~XStimParams.HiDynamicRange
    BN = questdlg('Use HighDynamic Range (preferable)?','Dynamic Range is set to low','High','Low (normal)','High');
    if ~isempty(findstr(BN,'Low'))    scaleFactor = TDT.scaleFactor;      end
end
nEpochmsecs = XStimParams.epoch_duration(1);       % durations for epochs1 and 2

% increment testnumber
if(exist1('H.space4fig') & get(H.space4_recorddata,'Value'))
    update_dataFN;
end
fclose all;

XStimParams.curr_stimdur = str2num(get(H.space4_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space4_ISI,'String'));
XStimParams.numreps = str2num(get(H.space4_numreps,'String'));
XStimParams.reset_flag = 0;

% params for XCorr between freq bands
XC = .5;
XC_tol = .1;


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
str1 = 'fc';
disp('This is a FULLY CUED space4 test')

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
eval(['save ' FN.temp_stim_path 'XStimParams_space4_' str1 ' Temp_params;'])
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
if size(XStimParams.locations,2)~=1
    B = questdlg('Push OK, pick ONLY ONE location. Then <CR>.','Location not picked','OK','OK','OK');
    pause;
end
locind = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,1) &...
    hrtfdirmat(2,:) == XStimParams.locations(2,1)));

% update display
setInfo_space4;
if get(H.space4_recorddata,'Value')
    tempstr = '     Nagel&Doupe stim  epochs:(';
    for iepoch = 1:nEpochs
        tempstr = [tempstr  num2str(XStimParams.mu(iepoch)) '+/-' num2str(XStimParams.sigma(iepoch)) ];
    end
    tempstr = [tempstr  ') ' num2str(XStimParams.epoch_duration(1)) 'msec   atten: ' num2str(XStimParams.curr_ABL)];
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
set(H.space4_status,'BackgroundColor','yellow');
set(H.space4_status,'ForegroundColor','blue');
set(H.space4_status,'String','Status: building stimulus');


if strcmp(XStimParams.test_type,'Space4_limitXC')
    set(H.stim_type,'Value',8);
end
stimtype = get(H.stim_type,'Value');

clear ENV STIM
trial_next = [];


if stimtype ==9         % from file
    fid = fopen([FN.stim_path FN.stim],'r');
    trial_next = fread(fid,inf,'float');
    fclose(fid);
    trial_next = trial_next - mean(trial_next);
    %trial_next = trial_next/mom(trial_next,2);
    
    % note scale factor changes with XStimParams.HiDynamicRange
    trial_next = trial_next(:)' * scaleFactor;
    
    %Add in the leading silent, trailing silent, and HRTFpts period
    trial_next =  [zeros(1,silence_lead_pts) trial_next zeros(1,silence_trail_pts + 254)];
  
    %%%%%%%%%%%%%%
elseif   strcmp(XStimParams.test_type,'Space4_limitXC')
    cF = 1000:100:11000;
    n_cF = length(cF);
    Fenv = 100;
    BBnoi = flatten(MakeBBNoise(TDT.Fs,nEpochmsecs+20));         % pad to eliminate ramping in ERBfilt
    E1 = zeros(n_cF,nEpochmsecs*30);
    for iepoch = 1:nEpochs
        eval(['SS' num2str(iepoch) ' = zeros(1,nEpochmsecs*30);']);
    end
    
    factor = .5;
    for icF = 1:n_cF
        NBnoi = use1_ERBfilt(BBnoi,TDT.Fs,cF(icF));
        NBnoi = NBnoi(301:end-300);
        NBnoi = NBnoi - mean(NBnoi);
        NBnoi = NBnoi / std(NBnoi);
        
        env0 = m_noi_exp(Fenv,nEpochmsecs/1000,TDT.Fs);
        env0 = norm2unif(env0);
        
        if icF > 1
            XCdiff = 1;
            while XCdiff > XC_tol
                env1 = E1(icF-1,:)*(1-factor) + env0*factor;
                XCcoef = corrcoef(E1(icF-1,:),env1);
                XCdiff = XC - XCcoef(1,2);
                if XCdiff>0
                    factor = factor - .01;
                else
                    factor = factor + .01;
                end
                XCdiff = abs(XCdiff);
            end
        else
            env1 = env0;
        end
        E1(icF,:) = env1;
        
        for iepoch = 1:nEpochs
            %temp = 10 .^ ((XStimParams.sigma(iepoch) .* env1 - XStimParams.sigma(iepoch)) /20);
            temp =  10 .^((XStimParams.mu(iepoch) + XStimParams.sigma(iepoch) .*env1)/20)-1;
            ind = find(temp<0);
            temp(ind) = zeros(size(ind));
            temp = NBnoi .* temp;
            eval(['SS' num2str(iepoch) ' = SS' num2str(iepoch) ' + temp;'])
        end
    end
    % renormalize to power==1
    for iepoch = 1:nEpochs
        eval(['SS' num2str(iepoch) ' = SS' num2str(iepoch) '/ std(SS' num2str(iepoch) ');'])
        eval(['trial_next = [trial_next SS' num2str(iepoch) '];']);
    end
    
    % ramp
    trial_next = ramp_sound(trial_next,TDT.Fs,5);
    
    %Add in the leading silent, trailing silent, and HRTFpts period
    trial_next =  [zeros(1,silence_lead_pts) trial_next zeros(1,silence_trail_pts + 254)];
    
    if get(H.space4_recorddata,'Value')
        % save sound to disk (for off-line STA analysis)
        eval(['save ' FN.temp_stim_path2 FN.data '_REP_' num2str(irep+1) '.space4noi trial_next']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    
    
else
    % make each epoch's carrier and envelope
    for iepoch = 1:nEpochs
        % carrier
        switch stimtype    
            case 1                  % tone
                stim = MakeTone(TDT.Fs,XStimParams.curr_freq,nEpochmsecs,0);
            case 2                  % gammatone
                stim = MakeBBNoise(TDT.Fs,nEpochmsecs);
                stim = use1_ERBfilt(stim,TDT.Fs,XStimParams.curr_freq);
            case 6                %Narrowband Noise
                stim = MakeNBNoise(TDT.Fs,nEpochmsecs,XStimParams.bandlimits);
            case 8                %Broadband Noise
                stim = MakeBBNoise(TDT.Fs,nEpochmsecs);
            case 10
                stim = MakeBBNoise(TDT.Fs,nEpochmsecs);
                stim = flatten (stim);        
            otherwise
                set(H.stim_type,'Value',10);
                stim = MakeBBNoise(TDT.Fs,nEpochmsecs);
                stim = flatten (stim);        
                disp('Unsupported stimulus type for space test; reset to BROADBAND-LNN');
        end
        
        % remove DCoffset & normalize ACPower==1
        stim = stim - mom(stim,1);
        stim = stim / mom(stim,2);
        %scale stim
        stim = stim * scaleFactor;
        STIM(iepoch,:) = stim;
        % make envelope
        env = m_noi_exp(100,nEpochmsecs/1000,TDT.Fs);
        env = norm2unif(env);         % makes distribution of amplitues almost uniform
        env = env /max1(abs(env));
        env = env/2+.5;
        env =  10 .^((XStimParams.mu(iepoch) + XStimParams.sigma(iepoch) .*env)/20)-1;
        ENV(iepoch,:) = env;
        
        % combine envelope and carrier
        stim =  stim(:) .* env(:);
        
        % re-normalize to ACpower
        factor = 10^(XStimParams.mu(iepoch)/20) *scaleFactor / mom(stim,2);
        stim = stim * factor;
        
        trial_next = [trial_next stim(:)'];
    end     % iepoch
    
    % ramp
    if get(H.stim_type,'Value') ~= 9
        trial_next = ramp_sound(trial_next,TDT.Fs,5);
    end
    
    %Add in the leading silent, trailing silent, and HRTFpts period
    trial_next =  [zeros(1,silence_lead_pts) trial_next zeros(1,silence_trail_pts + 254)];
    
    if get(H.space4_recorddata,'Value')
        % save sound to disk (for off-line STA analysis)
        eval(['save ' FN.temp_stim_path2 FN.data '_REP_1.space4noi trial_next']);
    end
end     % choose stimtype==9 or other

next_buf = BUF.L1;
S232('push16',trial_next(:)',Npts_totalplay);
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
finalspikematrix = 0;
%loop for reps
set(H.space4_status,'String','Status: playing stimuli');
set(H.space4_status,'BackgroundColor','green');
set(H.space4_status,'ForegroundColor','white');
set(H.space4_remtrials,'String',num2str(0));

param3 = 1;     % keeps track of envelope used (either 1 or irep)
for irep = 1:XStimParams.numreps
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    irep_trial = [];
    Nspikes = [];
    
    %Check for pause by user
    if pause_check    return; end
    
    % load HRTFs
    %Load left      
    S232('pushf',hrTF1(locind,:),255);
    S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
    %Load right
    S232('pushf',hrTF2(locind,:),255);
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
    if irep < XStimParams.numreps
        trial_next = [];
        if stimtype ==9         % from file
            if rand(1,1) > XStimParams.focalProb
                param3(irep+1) = round(rand(1,1)*9+1.4999999);
            else
                param3(irep+1) = 1;
            end
            ind = findstr(FN.stim,'.');
            tempFN = [FN.stim(1:ind-2) num2str(param3(irep+1)) FN.stim(ind:end)]
            
            fid = fopen([FN.stim_path tempFN],'r');
            trial_next = fread(fid,inf,'float');
            fclose(fid);
            trial_next = trial_next - mean(trial_next);
            trial_next = trial_next/mom(trial_next,2);
            % note scale factor changes with XStimParams.HiDynamicRange
            trial_next = trial_next(:)' * scaleFactor;
            
            %Add in the leading silent, trailing silent, and HRTFpts period
            trial_next =  [zeros(1,silence_lead_pts) trial_next zeros(1,silence_trail_pts + 254)];
            
            
         %%%%%%%%%%%% start of new   
        elseif   strcmp(XStimParams.test_type,'Space4_limitXC')
            E1 = zeros(n_cF,nEpochmsecs*30);
            for iepoch = 1:nEpochs
                eval(['SS' num2str(iepoch) ' = zeros(1,nEpochmsecs*30);']);
            end
            
            factor = .5;
            for icF = 1:n_cF
                NBnoi = use1_ERBfilt(BBnoi,TDT.Fs,cF(icF));
                NBnoi = NBnoi(301:end-300);
                NBnoi = NBnoi - mean(NBnoi);
                NBnoi = NBnoi / std(NBnoi);
                
                env0 = m_noi_exp(Fenv,nEpochmsecs/1000,TDT.Fs);
                env0 = norm2unif(env0);
                
                if icF > 1
                    XCdiff = 1;
                    while XCdiff > XC_tol
                        env1 = E1(icF-1,:)*(1-factor) + env0*factor;
                        factor = factor+.01;
                        XCcoef = corrcoef(E1(icF-1,:),env1);
                        XCdiff = XC - XCcoef(1,2);
                        if XCdiff>0
                            factor = factor - .01;
                        else
                            factor = factor + .01;
                        end
                        XCdiff = abs(XCdiff);
                    end
                else
                    env1 = env0;
                end
                E1(icF,:) = env1;
                
                for iepoch = 1:nEpochs
                    temp =  10 .^((XStimParams.mu(iepoch) + XStimParams.sigma(iepoch) .*env1)/20)-1;
                    ind = find(temp<0);
                    temp(ind) = zeros(size(ind));
                    temp = NBnoi .* temp;
                    eval(['SS' num2str(iepoch) ' = SS' num2str(iepoch) ' + temp;'])
                end
            end
            % renormalize to power==1
            for iepoch = 1:nEpochs
                eval(['SS' num2str(iepoch) ' = SS' num2str(iepoch) '/ std(SS' num2str(iepoch) ');'])
                eval(['trial_next = [trial_next SS' num2str(iepoch) '];']);
            end
            
            % ramp
            trial_next = ramp_sound(trial_next,TDT.Fs,5);
            
            %Add in the leading silent, trailing silent, and HRTFpts period
            trial_next =  [zeros(1,silence_lead_pts) trial_next zeros(1,silence_trail_pts + 254)];
            
            if get(H.space4_recorddata,'Value')
                % save sound to disk (for off-line STA analysis)
                eval(['save ' FN.temp_stim_path2 FN.data '_REP_' num2str(irep+1) '.space4noi trial_next']);
            end
 
         %%%%%%%%%%%%%%%%%% end of new   
            
            
        else
            
            if rand(1,1) > XStimParams.focalProb
                param3(irep+1) = irep+1;
            else
                param3(irep+1) = 1;
            end
            for iepoch = 1:nEpochs
                if param3(irep+1) >1
                    switch stimtype    
                        case 1                  % tone
                            stim = MakeTone(TDT.Fs,XStimParams.curr_freq,nEpochmsecs,0);
                        case 2                  % gammatone
                            stim = MakeBBNoise(TDT.Fs,nEpochmsecs);
                            stim = use1_ERBfilt(stim,TDT.Fs,XStimParams.curr_freq);
                        case 6                %Narrowband Noise
                            stim = MakeNBNoise(TDT.Fs,nEpochmsecs,XStimParams.bandlimits);
                        case 8                %Broadband Noise
                            stim = MakeBBNoise(TDT.Fs,nEpochmsecs);
                        case 10
                            stim = flatten(MakeBBNoise(TDT.Fs,nEpochmsecs));
                    end
                    % remove DCoffset & normalize ACPower==1
                    stim = stim - mom(stim,1);
                    %stim = stim / mom(stim,2);
                    %scale stim
                    stim = stim * scaleFactor;
                    % make envelope
                    env = m_noi_exp(100,nEpochmsecs/1000,TDT.Fs);
                    env = norm2unif(env);         % makes distribution of amplitues almost uniform
                    env = env /max1(abs(env));
                    env = env/2+.5;
                    env =  10 .^((XStimParams.mu(iepoch) + XStimParams.sigma(iepoch) .*env)/20)-1;
                else
                    env = ENV(iepoch,:);
                    stim = STIM(iepoch,:);    
                end
                % combine envelope and carrier
                stim =  stim(:) .* env(:);
                
                % re-normalize to ACpower
                factor = 10^(XStimParams.mu(iepoch)/20) *scaleFactor / mom(stim,2);
                stim = stim * factor;
                
                trial_next = [trial_next stim(:)'];
            end     % iepoch
            
            % ramp
            if get(H.stim_type,'Value') ~= 9
                trial_next = ramp_sound(trial_next,TDT.Fs,5);
            end
            
            %Add in the leading silent, trailing silent, and HRTFpts period
            trial_next =  [zeros(1,silence_lead_pts) trial_next zeros(1,silence_trail_pts + 254)];
            
            if get(H.space4_recorddata,'Value')
                % save sound to disk (for off-line STA analysis)
                eval(['save ' FN.temp_stim_path2 FN.data '_REP_' num2str(irep+1) '.space4noi trial_next']);
            end
        end             % stimval = 9?
        
        % pop sound to DAMA
        if next_buf == BUF.L1   
            next_buf = BUF.R1;  next_seq = BUF.playseq_R;
        else
            next_buf = BUF.L1;  next_seq = BUF.playseq_L;
        end
        S232('push16',trial_next,Npts_totalplay);
        S232('qpop16',next_buf);     
    end             % repnum-1?
    
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
    if(exist1('H.space4fig') & ~isempty(spikes)) 
        spikes = spikes(ind);
        spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
        EL_trial = [EL_trial;hrtfdirmat(1,locind) * ones(size(spikes))];
        AZ_trial = [AZ_trial;hrtfdirmat(2,locind) * ones(size(spikes))];
        irep_trial = [irep_trial;irep * ones(size(spikes))];
        Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
    end
    finalspikematrix = finalspikematrix + length(spikes);
    if pause_check    return; end
    
    %%%%%% reMake play specification buffer
    S232('dpush',10);
    S232('value',0);
    S232('make',0,next_seq);
    S232('make',1,next_seq);
    S232('make',2,0);
    S232('qpop16',BUF.playspec);
    %%%%%%%%%
    
    while t < ISI/1000
        t = toc;
        if pause_check    return; end
    end
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/irep;
    if(exist1('H.space4fig') & ~exist1('H.space4_finalspikeratefig'))
        H.space4_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','Space Test Spike Rate Plot',...
            'NumberTitle','off');
        H.space4_spikeaxes = axes;
    end
    figure(H.space4_finalspikeratefig)
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.space4_spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(irep)]);
    colorbar
    
    %Record Data
    datamatrix = [datamatrix;[Nspikes spikes_trial irep_trial EL_trial AZ_trial]];
    if(exist1('H.space4fig') & get(H.space4_recorddata,'Value'))
        record_data3(XStimParams,datamatrix,param3);
    end
    
    set(H.space4_remreps,'String',num2str(XStimParams.numreps - irep));
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.space4_finalspikeratefig)
set(H.space4_finalspikeratefig,'Name','Final Plot for Space Test');
plotdiam1(XStimParams.locations, interimspikerate);
set(H.space4_spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
colorbar

% plot all epochs
% still to do ????

set(H.space4_status,'String','Status: results');
set(H.space4_status,'BackgroundColor','blue');

set(H.exitspace4,'Visible','on');
set(H.resetspace4,'Visible','on');

% increment test number
if(exist1('H.space4fig') & get(H.space4_recorddata,'Value'))
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
                    set(H.space4_modfile,'String',[FN.mod_path FN.mod]);
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
while (exist1('H.space4fig') & get(H.pausespace4,'Value'))
    pause(0);
    if(~exist1('H.space4fig')) return; end         
    set(H.exitspace4,'Visible','on');
    set(H.resetspace4,'Visible','on');
    if(exist1('H.space4fig') & get(H.resetspace4,'Value') == 1)
        set(H.resetspace4,'Value',0);
        set(H.pausespace4,'Value',0);
        Reset_space;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_space4;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end