function [] = Engage_Delay3()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_Delay

%*******************************************************************************
% Delay 
%*******************************************************************************
% NOTE: (June 2, 2005) Copied from Kip's 2Source
% revising using PStim November 10, 2005

%Put parameters into XStimParams
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.reset_flag = 0;

[TrialSeq, ColumnLabels, DelayModLabels] = GetDelayTrialSeq;
%assignin('base','TrialSeq',TrialSeq);

PStimG = DelayParams2PStimG(XStimParams, FN); % save PStimG instead of XStimParams and FN

Trials = size(TrialSeq,1);
MaxVars = size(TrialSeq,2);
% make data matrix and add a row (or 2) for NumSpikes
RLdata = [TrialSeq zeros(size(TrialSeq,1), 2)];
repNum = 0; % count rep 1 on initial NaN in trial sequence


%check if filt files assigned
if exist([FN.ephone_path FN.ephone2]) ~= 2
    ephonefilediagbox;
end
if XStimParams.space_flag
    while exist1([FN.space_path FN.space]) ~=2
        [FN.space,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select Fully-cued HRTF File');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space]);
        end
        FN.HRTFfiletype(1) = testHRTFfiletype(FN.space_path, FN.space);
    end
else
    if XStimParams.space_flag == 0
        disp('No HRTF filtering')
    end
end
		
%Specify DAMA buffers
clear BUF
BUF.L1				= 1;
BUF.R1				= 2;
BUF.playseq_L1		= 3;
BUF.playseq_R1		= 4;
BUF.playspec1		= 5;

%Make play sequence buffers
S232('allot16',BUF.playseq_L1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.L1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_L1);

S232('allot16',BUF.playseq_R1,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.R1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_R1);

%Make play specification buffer
S232('allot16',BUF.playspec1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L1);
S232('make',1,BUF.playseq_R1);
S232('make',2,0);
S232('qpop16',BUF.playspec1);

SR = TDT.Fs;
        
%Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
silence_len = (XStimParams.silence_lead * round(SR/1000));
%Add a piece of silence after stimulus 
silence_len2 = (XStimParams.silence_trail * round(SR/1000));

%Make Stimulus buffers
DUR = XStimParams.curr_stimdur;

% assuming delay + lag segment extension will not exceed DUR (hence 2*DUR)
TotalStimPnts = silence_len2 + silence_len + 2*DUR*round(SR/1000) + TDT.ephonefiltlen + TDT.hrtffiltlen;
TotalStimDur = (TotalStimPnts/SR)*1000;
S232('allot16',BUF.L1, TotalStimPnts);
S232('allot16',BUF.R1, TotalStimPnts);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/SR);
S232('PD1npts',1, TotalStimPnts);

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


%Set MII parameters
mii_us_per_sample = 10;     %microsecond per sample
mii_separation = 100; 		%only take events separated by 100 samples (i.e., 1 ms)

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

ISI = XStimParams.test_ISI;

RasterRows = 30;
%RasterData = zeros(RasterRows, TrialInRep?);

% get HRTFs before generating signals as this should not change        
% get HRTFs from variables (could get from sequence...)

dir = 0;
if FN.HRTFfiletype(1) == 2
    eval(['load -mat ' FN.space_path FN.space]);
    TF1_space = TF1; TF2_space = TF2;
    dir_space = dir;
end
if FN.HRTFfiletype(1) == 1
    hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space]));
elseif FN.HRTFfiletype(1) == 2
    hrtfdirmat = dir_space;
else
    warning('unknown HRTF file type'); return;
end
locind_1 = max(find(hrtfdirmat(1,:) == XStimParams.locations(2) & hrtfdirmat(2,:) == XStimParams.locations(1)));
if FN.HRTFfiletype(1) == 1
    Lref_1 = mtlrch([FN.space_path FN.space],(2*locind_1)-1);
    Rref_1 = mtlrch([FN.space_path FN.space],2*locind_1);
else
    Lref_1 = TF1_space(locind_1,:)';
    Rref_1 = TF2_space(locind_1,:)';
end
locind_2 = max(find(hrtfdirmat(1,:) == XStimParams.offset_el(2) & hrtfdirmat(2,:) == XStimParams.offset_az(2)));
if FN.HRTFfiletype(1) == 1
    Lref_2 = mtlrch([FN.space_path FN.space],(2*locind_2)-1);
    Rref_2 = mtlrch([FN.space_path FN.space],2*locind_2);
else
    Lref_2 = TF1_space(locind_2,:)';
    Rref_2 = TF2_space(locind_2,:)';
end


TrialCnt=1;
RowNum=1;

set(H.status,'String','Status: Playing Stimuli');
set(H.status,'BackgroundColor','green');
set(H.remreps,'String',num2str(XStimParams.numreps));

% increment testnumber  ???
if(exist1('H.Delayfig') & get(H.recorddata,'Value'))
    update_dataFN;
end

% loop for trials
tic
while (exist1('H.Delayfig') & (TrialCnt <= Trials))
    
    %Check for pause by user
    if pause_check  return; end
    
    if isnan(TrialSeq(TrialCnt,3))
        %pause(1); % pause to allow for adjustments...
        repNum=repNum+1; % counts reps when NaN found
    end

    % ????
    set(H.remreps,'String',num2str(XStimParams.numreps-repNum));
    TrialsRep = Trials/XStimParams.numreps;
    set(H.remtrials,'String',num2str(TrialsRep-(TrialsRep-TrialCnt)-(TrialsRep*(repNum-1))));
    set(H.remtrials,'String',num2str((TrialsRep*repNum)-TrialCnt));
    %set(DelayCtrl.TrialRemaining, 'String', [num2str(TrialsRep-(TrialsRep-TrialCnt)-(TrialsRep*(repNum-1))) ' of ' num2str(TrialsRep)]);
   
    tDel = TrialSeq(TrialCnt,9);
    
    SingleStim = 0;
    if tDel > DUR*1.5
        tDel = 0;
        SingleStim = 2; % stim #2 alone
    elseif tDel < -DUR*1.5
        tDel = 0;
        SingleStim = 1; % stim #1 alone
    end
    

    if ~isnan(TrialSeq(TrialCnt,3)) % skip if not playing a sound
        OnsetDelay_1=0; OnsetDelay_2=0; OffsetDelay_1=0; OffsetDelay_2=0;
        if TrialSeq(TrialCnt,13) == 3 | TrialSeq(TrialCnt,13) == 5 % no lead segment
            if tDel > 0     % Stim#1 lags
                OnsetDelay_2 = -abs(tDel);
            else            % Stim#2 lags
                OnsetDelay_1 = -abs(tDel) + TrialSeq(TrialCnt,12);
            end
	
        elseif TrialSeq(TrialCnt,13) == 4 | TrialSeq(TrialCnt,13) == 6 % no lag segment
            if tDel > 0     % Stim#1 lags
                OffsetDelay_1 = -abs(tDel);
            else            % Stim#2 lags
                OffsetDelay_2 = -abs(tDel) + TrialSeq(TrialCnt,12);
            end
        else
            if tDel > 0     % Stim#1 lags
                OffsetDelay_1 = TrialSeq(TrialCnt,12);
            else            % Stim#2 lags
                OffsetDelay_2 = TrialSeq(TrialCnt,12);
            end
        end
        
        OnsetPnts_1=round((OnsetDelay_1/1000)*SR); % Time to Onset
        OnsetPnts_2=round((OnsetDelay_2/1000)*SR);
        OffsetPnts_1=round((OffsetDelay_1/1000)*SR); % Time to Offset
        OffsetPnts_2=round((OffsetDelay_2/1000)*SR);
	
        % duration of stimulus in points adjusted for Onset and Offset envelope shift/delay times
        Ch1_Pnts = round((DUR/1000)*SR) + OnsetPnts_1 + OffsetPnts_1;
        Ch2_Pnts = round((DUR/1000)*SR) + OnsetPnts_2 + OffsetPnts_2;
        
        % onset and offset ramp duations
        whichRamp=get(H.DelayOnOff_mode,'Value'); % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
        if whichRamp == 1 & tDel >= 0 % stim #1 lags
            OnsetMs_1 = TrialSeq(TrialCnt,10);
            OffsetMs_1 = TrialSeq(TrialCnt,11);
            OnsetMs_2 = XStimParams.ramp_timeS_on(1,1);
            OffsetMs_2 = XStimParams.ramp_timeS_off(1,1);
        elseif whichRamp == 1 & tDel < 0 % stim #1 leads
            OnsetMs_1 = XStimParams.ramp_timeS_on(1,1);
            OffsetMs_1 = XStimParams.ramp_timeS_off(1,1);
            OnsetMs_2 = TrialSeq(TrialCnt,10);
            OffsetMs_2 = TrialSeq(TrialCnt,11);
        elseif whichRamp == 2
            OnsetMs_1 = TrialSeq(TrialCnt,10);
            OffsetMs_1 = TrialSeq(TrialCnt,11);
            OnsetMs_2 = TrialSeq(TrialCnt,10);
            OffsetMs_2 = TrialSeq(TrialCnt,11);
        elseif whichRamp == 3  & tDel >= 0 % stim #1 lags
            OnsetMs_1 = XStimParams.ramp_timeS_on(1,1);
            OffsetMs_1 = XStimParams.ramp_timeS_off(1,1);
            OnsetMs_2 = TrialSeq(TrialCnt,10);
            OffsetMs_2 = TrialSeq(TrialCnt,11);
        elseif whichRamp == 3  & tDel < 0 % stim #1 leads
            OnsetMs_1 = TrialSeq(TrialCnt,10);
            OffsetMs_1 = TrialSeq(TrialCnt,11);
            OnsetMs_2 = XStimParams.ramp_timeS_on(1,1);
            OffsetMs_2 = XStimParams.ramp_timeS_off(1,1);
        else
            warning('No stimulus ramp!');
        end
        % set above 2000*(1/SR) or 2 samples to avoid dividing by zero
        OnsetMs_1 = max(OnsetMs_1, 2000*(1/SR));
        OffsetMs_1 = max(OffsetMs_1, 2000*(1/SR));
        OnsetMs_2 = max(OnsetMs_2, 2000*(1/SR));
        OffsetMs_2 = max(OffsetMs_2, 2000*(1/SR));
        
        % Stimulus envelopes
        Stim1_env = ones(1,Ch1_Pnts);
        Stim1_env = ramp_sound_onoff(Stim1_env, SR, OnsetMs_1, OffsetMs_1);
        Stim2_env = ones(1,Ch2_Pnts);
        Stim2_env = ramp_sound_onoff(Stim2_env,SR, OnsetMs_2, OffsetMs_2);
        
        % reinsert onset envelope points
        if TrialSeq(TrialCnt,13) == 3 | TrialSeq(TrialCnt,13) == 5 % no lead segment
            if tDel > 0     % Stim#1 lags
                Stim2_env = [zeros(1,abs(OnsetPnts_2)) Stim2_env];
            else                            % Stim#2 lags
                Stim1_env = [zeros(1,abs(OnsetPnts_1)) Stim1_env];
            end
        end
                  
        % make long (300 ms) noise
        Stim1_L = GetNoise_BB(TrialSeq(TrialCnt,7), SR, 300);
	
        % if uncorrelated
        if TrialSeq(TrialCnt,13) == 2 | TrialSeq(TrialCnt,13) == 5 | TrialSeq(TrialCnt,13) == 6
            Stim2_L = fliplr(Stim1_L);
        else
            Stim2_L = Stim1_L;
        end
        
        % crop to duration of envelope
        Stim1_L = Stim1_L(1,1:size(Stim1_env,2));
        Stim2_L = Stim2_L(1,1:size(Stim2_env,2));
        
        % length (points) before convolution
        Ch1_Pnts = size(Stim1_L,2);
        Ch2_Pnts = size(Stim2_L,2);
        
        % remove any DCoffset
        Stim1_L = Stim1_L - mom(Stim1_L,1);
        Stim2_L = Stim2_L - mom(Stim2_L,1);
        
        % Envelope stimuli
        Stim1_L = Stim1_L .* Stim1_env;
        Stim2_L = Stim2_L .* Stim2_env;
         
        % silence one channel if playing only one stimulus
        % (i.e., delay exceeds stimulus duration * 1.5)
        if SingleStim == 1
            Stim2_L = Stim2_L*0;
        elseif SingleStim == 2
            Stim1_L = Stim1_L*0;
        end
	
        % convolve stim #1 with HRTFs
        Stim1_R=Stim1_L;
        Stim1_L = conv(Stim1_L,Lref_1);
        Stim1_R = conv(Stim1_R,Rref_1);
        
        % convolve stim #2 with HRTFs
        Stim2_R=Stim2_L;
        Stim2_L = conv(Stim2_L,Lref_2);
        Stim2_R = conv(Stim2_R,Rref_2);
                    
        % Delay by inserting points
        if tDel > 0     % Stim#1 lags
            Stim1_L = [zeros(1,round((abs(tDel)/1000)*SR)) Stim1_L];
            Stim1_R = [zeros(1,round((abs(tDel)/1000)*SR)) Stim1_R];
        else            % Stim#2 lags
            Stim2_L = [zeros(1,round((abs(tDel)/1000)*SR)) Stim2_L];
            Stim2_R = [zeros(1,round((abs(tDel)/1000)*SR)) Stim2_R];
        end
        
        % make Ch1 stimuli have equal points = TotalStimPnts
		Stim1_L = [Stim1_L zeros(1,TotalStimPnts-size(Stim1_L,2))];
		Stim2_L = [Stim2_L zeros(1,TotalStimPnts-size(Stim2_L,2))];
       
        Stim1_R = [Stim1_R zeros(1,TotalStimPnts-size(Stim1_R,2))];
        Stim2_R = [Stim2_R zeros(1,TotalStimPnts-size(Stim2_R,2))];
        
        % plot each left and right stimuli?
        if 0 %& TrialCnt == TrialToPlot
            Stim1_Delay =0; Stim2_Delay =0; 
            if tDel > 0 % Stim#1 lags
                Stim1_Delay = abs(tDel);
            else        % Stim#2 lags
                Stim2_Delay = abs(tDel);
            end
            PlotDelayedStim(Stim1_env, Stim2_env, Stim1_L, Stim1_R, Stim2_L, Stim2_R, Stim1_Delay, Stim2_Delay);
        end
        
        % Combine sounds
        Ch1_Snd = (Stim1_L + Stim2_L)/2;
        Ch2_Snd = (Stim1_R + Stim2_R)/2;
        
        % re-scale
        ABAval = 0.5*(mom(Ch1_Snd,2) + mom(Ch2_Snd,2));
        scalefact = TDT.scalevalue/ABAval;
        Ch1_Snd = round(scalefact*Ch1_Snd);
        Ch2_Snd = round(scalefact*Ch2_Snd);
        
    else % make zeros
        Ch1_Snd = GetNoise_BB(1, SR, TotalStimDur)*0;
        Ch2_Snd = Ch1_Snd;
    end
    

    %Wait till PD1 is finished
    while S232('PD1status',1) usec_delay(1000); end
    
    % Load stimuli
    S232('push16',Ch1_Snd,length(Ch1_Snd));
    S232('qpop16',BUF.L1);
    S232('push16',Ch2_Snd,length(Ch2_Snd));
    S232('qpop16',BUF.R1);
    
    %Check for pause by user
    if pause_check  return; end

    if(exist1('H.Delayfig') & (TrialCnt <= Trials))           
        S232('seqplay',BUF.playspec1);
        S232('PD1arm',1);
        
        %Send trigger
        %Set up MII
        m100x( C_.INIT );
        m110dx( C_.INIT );
        m110dx( C_.CLOCK, mii_us_per_sample);
        m110dx( C_.MODE, M_.PST );
        
        if (TrialCnt <= Trials)
            while toc < ISI/1000     end
            %Start clock
            m110dx( C_.START);
            %Send pulse: PD1 GO!
            m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
            tic
        end
        
    end

    %Wait until PD1 is finished
    while S232('PD1status',1) usec_delay(1000); end
    S232('PD1stop',1);
    
    %Stop the m110 and get spikes
    if(TrialCnt >= 1)							% first trial NOT just for loading
        m110dx( C_.STOP);
        spikes = m110dx( C_.DATA, 3*DUR); 			% Take XStimParams.curr_stimdur spikes max +50 ms!!!!
        ind = find(spikes ~= 0); 						% Get clock events that are spikes
        spikes = spikes(ind);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
        spikes = spikes(ind);
        spikes = spikes/(1000/mii_us_per_sample);
        spikes = spikes';
        
        if(exist1('H.Delayfig') & ~isempty(spikes)) 
            if RowNum > RasterRows
                RowNum=1;
            end
            
            if 0 %PStimG.SimpleRaster % plot spikes
                %disp([RasterRows RowNum spikes]);
                PStim_SRaster(PStimG, RasterRows,RowNum, spikes);
            end
            
            numSpikes = length(spikes);
             if size(RLdata,2) < MaxVars+2+numSpikes
                 RLdata = [RLdata zeros(size(RLdata,1), (MaxVars+2+numSpikes)-size(RLdata,2))];
             end
           
            RLdata(TrialCnt, MaxVars+1) = numSpikes;
            RLdata(TrialCnt, MaxVars+2 : MaxVars+2+numSpikes-1) = spikes;
            
        end
    end
        

    %Record Data
    if exist1('H.Delayfig') & TrialsRep-(TrialsRep-TrialCnt)-(TrialsRep*(repNum-1)) == TrialsRep % last trial
    
        DelayData = RLdata(1:TrialCnt,:);
        if get(H.recorddata,'Value')
            eval(['save ' [FN.data_path FN.data '.mat']  ' DelayData ColumnLabels PStimG;']) % save locations
        end
        PStim_TestRaster(DelayData(1:TrialCnt,:), [13 9], PStimG);
        Delay_PlotData(DelayData(1:TrialCnt,:), PStimG);
    end
    
    TrialCnt=TrialCnt+1;
    RowNum=RowNum+1;
    
    pause(0);
end %end loop over trials


set(H.status,'BackgroundColor','blue');
set(H.status,'String','Status: Results');
set(H.exitDelay,'Visible','on');
set(H.resetDelay,'Visible','on');

% increment test number
if(exist1('H.Delayfig') & get(H.recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end


% //////////////////////////////////////////////////
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
	while (exist1('H.Delayfig') & get(H.pauseDelay,'Value'))
        pause(0);
        if(~exist1('H.Delayfig')) return; end         
        set(H.exitDelay,'Visible','on');
        set(H.resetDelay,'Visible','on');
        if(exist1('H.Delayfig') & get(H.resetDelay,'Value') == 1)
            set(H.resetDelay,'Value',0);
            set(H.pauseDelay,'Value',0);
            Reset_Delay;   flag = 1;
            return;
        end
        if isempty(XStimParams.locations)
            Reset_Delay;   flag=1;
            return;
        end
	end
	
	if XStimParams.reset_flag
        flag = 1;
        XStimParams.reset_flag = 0;
	end

    
% //////////////////////////////////////////////////
function [bbnoise] = GetNoise_BB(state,Fs,dur)
    %GetNoise_BB:	Create a BroadBand Noise (2-11 kHz)
    %state:  rand state
    %Fs:	Sampling rate (Hz)
    %dur:	Stimulus duration (ms)
    dur = dur/1000;
    len = dur*Fs;
    minfreq = 2000; maxfreq = 11000;
    minfreq = round(((minfreq + 1)/Fs) * len);
    maxfreq = round(((maxfreq + 1)/Fs) * len);
    range = maxfreq-minfreq+1;
    % mag spectrum = 1 between set frequencies:
    mag = zeros(len,1);
    mag(minfreq:maxfreq) = ones(range,1);
    % random phase spectrum between set frequencies:
    % rand('state',sum(100*clock));
    rand('state',state); % use state
    phi = (rand(len,1) - 0.5) * (2*pi);
    % combine phase and magnitude:
    X = mag .* ( (cos(phi)) + (i .* sin(phi)) );
    % convert to time domain:
    bbnoise = real(ifft(X));
    bbnoise = (bbnoise/max(abs(bbnoise)))';

         
% //////////////////////////////////////////////////
function [ramped_sound] = ramp_sound_onoff(var_name,Fs,ramp_time_on, ramp_time_off)
	%Ramp_Sound, [ramped_sound] = ramp_sound_onoff(var_name,Fs,ramp_time)
	%var_name:		sound to be ramped
	%Fs:			sampling rate (Hz)
	%ramp_time_on:		onset time in ms
	%ramp_time_off:		offset time in ms
	
	if(~isreal(var_name))
       error('Sound variable must contain a real array');
	else
       [r,c] = size(var_name);
       if(r > 1 & c > 1)
          error('Sound variable must be a 1-D arrary');
       end
	end
	
	lenvarname = length(var_name);
	ramp = ones(size(var_name));
	samp_per = 1000/Fs; %in ms
	num_pts_on = round(ramp_time_on/samp_per);
	num_pts_off = round(ramp_time_off/samp_per);
	if(num_pts_on>=3)
        ramp(1:num_pts_on) = 0:1/num_pts_on:1-(1/num_pts_on);
	end
	if(num_pts_off>=3)
        ramp(lenvarname - num_pts_off+1:lenvarname) = 1-(1/num_pts_off):-1/num_pts_off:0;
	end
	ramped_sound = var_name .* ramp;
	return

