function [] = Engage_LongSAMs()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_LongSAMs

%*******************************************************************************
%	The LongSAMs Test operation
% altered version to eliminate double buffering
% should use HRIRs without earphone equalization - usually called *.std
%*******************************************************************************

stimuli_dir = FN.temp_stim_path;
fclose all;

if XStimParams.HiDynamicRange
    scaleFactor = 60;
else
    scaleFactor = TDT.scaleFactor;
end

% reset stim_type to BroadBand
    set(H.stim_type,'Value',8);

XStimParams.curr_ABL = str2num(get(H.LongSAMs_ABL,'String'));
ABL = abs(XStimParams.curr_ABL);
XStimParams.curr_stimdur = str2num(get(H.LongSAMs_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.LongSAMs_ISI,'String'));
XStimParams.numreps = str2num(get(H.LongSAMs_numreps,'String'));
XStimParams.reset_flag = 0;

%Specify DAMA buffers
clear BUF
BUF.L1				= 1;
BUF.R1				= 2;
BUF.playseq_L1		= 5;
BUF.playseq_R1		= 6;
BUF.playspec1		= 7;

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

% Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
DUR_silence = XStimParams.silence_lead; 					%ms
silence_len = (DUR_silence * round(TDT.Fs/1000));
% Add a piece of silence after stimulus 
DUR_silence2 = XStimParams.silence_trail; 					%ms
silence_len2 = (DUR_silence2 * round(TDT.Fs/1000));


%Make Stimulus buffers
DUR = XStimParams.curr_stimdur;
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen + TDT.ephonefiltlen;
    XStimParams.test_type = 'LongSAMs FC';     str1 = 'fc';
    disp('This is a FULLY CUED LongSAMs test')

S232('allot16',BUF.L1,Npts_totalplay);
S232('allot16',BUF.R1,Npts_totalplay);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,Npts_totalplay);

%Get Earphone filters
if FN.HRTFfiletype(6) == 1
    ephonefname = [FN.ephone_path FN.ephone2];
    ephonefilt_left  = (mtlrch(ephonefname,1))';
    ephonefilt_right = (mtlrch(ephonefname,2))';
else
    eval(['load -mat ' FN.ephone_path FN.ephone2]);
    ephonefilt_left  = TF1;
    ephonefilt_right = TF2;
    clear TF1 TF2 dir
end

%Load Earphone filters
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

S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

ISI = XStimParams.test_ISI;
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000)); 				%correct for ITD filtlength
ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));			%correct for ephonefilt

%Get the HRTF indices
XStimParams.locations = GUI.locations1';
if ~length(XStimParams.locations)
    set(H.pickerfig,'Color', [.1 .2 .8]);
    set(H.picker_error,'visible','on');
    pause;
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end

% save XStimParams for this test
if get(H.LongSAMs_recorddata,'Value')
    tempstr = ['    ' str1 ' LongSAMs-test: ' num2str(size(XStimParams.locations,2)) ' locations;     atten: ' num2str(abs(ABL))];
    update_diary
end

Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_LongSAMs_' str1 ' Temp_params;'])
clear Temp_params str1

    if isempty(FN.space_std) | ~isempty(strfind(FN.space_std,'eq'))        % FN not yet picked
        [FN.space_std,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.std preferred)');
        if(FN.space_path ~= 0)
            set(H.LongSAMsfile,'String',[FN.space_path FN.space_std]);
        end
        set(H.LongSAMsfilepb,'Value',0);
        FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
    end
    %%%%%%%
    if FN.HRTFfiletype(1,2) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space_std]));
    elseif FN.HRTFfiletype(1,2) == 2
        dir = 0;
        eval(['load -mat ' FN.space_path FN.space_std]);
        TF1_LongSAMs = TF1;    TF2_LongSAMs = TF2;
        hrtfdirmat = dir;
        clear dir TF1 TF2
    else
        disp(['LongSAMs HRTFfiletype incorrect'])
        return
    end
clear locind
for locnum = 1:size(XStimParams.locations,2)
    locind(locnum) = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
end

% update display
setInfo_LongSAMs;


% update tempStruct for play2_record2b_SAMB
tempStruct.mdepth = XStimParams.mod_depth(1);
tempStruct.latten = abs(XStimParams.curr_ABL)-20;
tempStruct.ratten = abs(XStimParams.curr_ABL)-20;
tempStruct.SamplingPeriod = 1000000/TDT.Fs;
tempStruct.buf_pts = TDT.Fs / XStimParams.mod_freq(1);
tempStruct.nCycles_perRep = str2num(get(H.LongSAMstestmod_nCycles_perRep));
tempStruct.nptsTotalPlay = tempStruct.nCycles_perRep * tempStruct.buf_pts;
tempStruct.ISI = XStimParams.test_ISI;

%%%%% Loop to make the stimuli we'll use
remreps = 1;
set(H.LongSAMs_status,'String','Status: building stimuli');
set(H.LongSAMs_status,'BackgroundColor','blue');

set(H.LongSAMs_remreps,'String',num2str(remreps));
repnum = 1;
numtrials = size(XStimParams.locations,2);
finalspikematrix = zeros(1,numtrials);
clear stim_from_file

%Randomize the stimuli
randseq = randperm(numtrials);

trialnum = 1;
while (exist1('H.LongSAMsfig') & (trialnum <= numtrials))
    set(H.LongSAMs_status,'BackgroundColor','yellow');
    %Check for pause by user
    if pause_check    return; end
    
    %Make the stimuli
            %Broadband Noise
            [trial_left] = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
            trial_right = trial_left;
   
    % remove any DCoffset
    trial_left = trial_left - mom(trial_left,1);
    trial_right = trial_right - mom(trial_right,1);
    
    % normalize to ACPower (added 3/5/07 as in 2-source)
    if ~XStimParams.HiDynamicRange
        trial_left = trial_left / mom(trial_left,2);
        trial_right = trial_right / mom(trial_right,2);
    end

    
        %Ramp the stimuli
        ramp_time = 5; %ms
        [trial_left] = ramp_sound(trial_left,TDT.Fs,ramp_time);
        [trial_right] = ramp_sound(trial_right,TDT.Fs,ramp_time);
    
     % remove any DCoffset
    trial_left = trial_left - mom(trial_left,1);
    trial_right = trial_right - mom(trial_right,1);

    
    %Apply HRTF filtering
        if FN.HRTFfiletype(1,2) == 1
            hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locind(randseq(trialnum)))-1);
            hrtf_right = mtlrch([FN.space_path FN.space_std],2*locind(randseq(trialnum)));
        else
            hrtf_left = TF1_LongSAMs(locind(randseq(trialnum)),:);
            hrtf_right = TF2_LongSAMs(locind(randseq(trialnum)),:);
        end
        str1 = 'FC';
    
    trial_left = conv(trial_left, hrtf_left);
    trial_right = conv(trial_right, hrtf_right);

    % remove DC offset
    trial_left = trial_left - round(mean(trial_left));
    trial_right = trial_right - round(mean(trial_right));

    % adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
    % presented through DSPs (*.eq) at 0,0 (broadband)
    trial_left = trial_left * TDT.hrtf_Lfactor;
    trial_right = trial_right * TDT.hrtf_Rfactor;
        
    %scale stimuli 3/7/07
    % note scale factor changes with XStimParams.HiDynamicRange
    trial_left = trial_left * scaleFactor;
    trial_right = trial_right * scaleFactor;
    
    %Add in the leading silent period
    trial_left =  [zeros(1,silence_len) trial_left];
    trial_right = [zeros(1,silence_len) trial_right];
    
    %Add in the trailing silent period
    trial_left =  [trial_left zeros(1,silence_len2)];
    trial_right = [trial_right zeros(1,silence_len2)];
    
    %pad with zeros for ephonefilters
    trial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
    trial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
    
    % save stimuli to disk
    S232('push16',trial_left,length(trial_left));
    S232('qpop16',BUF.L1);
    fname = ['LongSAMs1_' str1 '.left_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
            '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
    evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
            [' ''' FN.temp_stim_path fname ''' '] ',0);'];
    eval(evalstr);
    %temp_left = dama2pc(BUF.L1);
    S232('push16',trial_right,length(trial_right));
    S232('qpop16',BUF.R1);
    fname = ['LongSAMs1_' str1 '.right_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
            '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
    evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
            [' ''' FN.temp_stim_path fname ''' '] ',0);'];
    eval(evalstr);
    %temp_right = dama2pc(BUF.R1);
    
    % set up for next trial                        
    remtrials = numtrials - trialnum;
    set(H.LongSAMs_remtrials,'String',num2str(remtrials));
    trialnum = trialnum + 1;
    set(H.LongSAMs_status,'BackgroundColor','blue');
    pause(0);
end 										%end loop over trials
%%%%%%%%%%%%%%%%%%% finished making stimuli
%Begin playing sounds   
set(H.LongSAMs_status,'String','Status: playing stimuli');
set(H.LongSAMs_status,'BackgroundColor','green');
set(H.LongSAMs_remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
buffcycle = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.LongSAMsfig') & get(H.LongSAMs_recorddata,'Value'))
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
%loop for reps
while (exist1('H.LongSAMsfig') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli
    randseq = randperm(numtrials);   
    trialnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    tic
    while (exist1('H.LongSAMsfig') & (trialnum <= numtrials) )
        %Check for pause by user
        if pause_check    return; end
        
        % load stimuli
        fname = ['LongSAMs1_' str1 '.left_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
                '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
        eval(evalstr);
        evalstr = ['S232(''qpop16'''  ',BUF.L1);'];
        eval(evalstr);
        fname = ['LongSAMs1_' str1 '.right_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
                '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
        eval(evalstr);
        evalstr = ['S232(''qpop16'''  ',BUF.R1);'];
        eval(evalstr);
                
        S232('seqplay',BUF.playspec1);
        S232('PD1arm',1);
        
        %Send trigger
        %Set up MII
        m100x( C_.INIT );
        m110dx( C_.INIT );
        m110dx( C_.CLOCK, mii_us_per_sample);
        m110dx( C_.MODE, M_.PST );
        
        while toc < ISI/1000 
            if pause_check    return; end
        end

        %Start clock
        m110dx( C_.START);
        %Send pulse: PD1 GO!
        m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
        tic
        
        pause(DUR/1000+.2);
        
        while(S232('PD1status',1)) usec_delay(1000);  end
        S232('PD1stop',1);
       
        %Stop the m110 and get spikes
        m110dx( C_.STOP);
        spikes = m110dx( C_.DATA, XStimParams.curr_stimdur*2); 			% Take 2*XStimParams.curr_stimdur spikes max
        ind = find(spikes ~= 0); 						% Get clock events that are spikes
        spikes = spikes(ind);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
        if(exist1('H.LongSAMsfig') & ~isempty(spikes)) 
            spikes = spikes(ind);
            spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
            EL_trial = [EL_trial;hrtfdirmat(1,locind(randseq(trialnum)))* ones(size(spikes))];
            AZ_trial = [AZ_trial;hrtfdirmat(2,locind(randseq(trialnum)))* ones(size(spikes))];
            repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
            Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
        end
        finalspikematrix(randseq(trialnum)) = finalspikematrix(randseq(trialnum)) + length(spikes);
        if pause_check    return; end
    
        remtrials = numtrials - trialnum;
        set(H.LongSAMs_remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end %end loop over trials
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.LongSAMsfig') & ~exist1('H.LongSAMs_finalspikeratefig'))
        H.LongSAMs_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','LongSAMs Test Spike Rate Plot',...
            'NumberTitle','off');
        H.LongSAMs_spikeaxes = axes;
    end
    figure(H.LongSAMs_finalspikeratefig)
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.LongSAMs_spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
    colorbar
    
    %Record Data
    if(exist1('H.LongSAMsfig') & get(H.LongSAMs_recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
        record_data3(XStimParams,datamatrix);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.LongSAMs_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.LongSAMs_finalspikeratefig)
set(H.LongSAMs_finalspikeratefig,'Name','Final Plot for LongSAMs Test');
plotdiam1(XStimParams.locations, interimspikerate);
set(H.LongSAMs_spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
        ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))], 'FontSize',8);
colorbar

set(H.LongSAMs_status,'String','Status: results');
set(H.LongSAMs_status,'BackgroundColor','blue');

set(H.exitLongSAMs,'Visible','on');
set(H.resetLongSAMs,'Visible','on');

% increment test number
if(exist1('H.LongSAMsfig') & get(H.LongSAMs_recorddata,'Value'))
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
                    set(H.LongSAMs_modfile,'String',[FN.mod_path FN.mod]);
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
while (exist1('H.LongSAMsfig') & get(H.pauseLongSAMs,'Value'))
    pause(0);
    if(~exist1('H.LongSAMsfig')) return; end         
    set(H.exitLongSAMs,'Visible','on');
    set(H.resetLongSAMs,'Visible','on');
    if(exist1('H.LongSAMsfig') & get(H.resetLongSAMs,'Value') == 1)
        set(H.resetLongSAMs,'Value',0);
        set(H.pauseLongSAMs,'Value',0);
        Reset_LongSAMs;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_LongSAMs;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end