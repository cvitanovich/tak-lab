function [] = Engage_space_RIR()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI
global RIR1
global RIR2

%Engage_space_RIR

%*******************************************************************************
%	The space_RIR Test operation
% no double buffering
% use HRIRs without earphone equalization - usually called *.std
% based on Engage_space_FPR on Oct 10,2012
%
%*******************************************************************************
colors = [ ...
    179 199 255; ...
    200 150 255; ...
    0    0  255; ...
    216 41  0; ...
    255 199 179;...
    255 150 200;...
    255   0    0;...
    199 255 179;...
    200 255 150;...
    0   255 0]/255;

fclose all;

if XStimParams.HiDynamicRange
    scaleFactor = 60;
else
    scaleFactor = TDT.scaleFactor;
end

XStimParams.curr_stimdur = str2num(get(H.space_RIR_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space_RIR_ISI,'String'));
XStimParams.numreps = str2num(get(H.space_RIR_numreps,'String'));
XStimParams.mod_freq(1) = str2num(get(H.space_RIR_mod_freq,'String'));
XStimParams.mod_depth(1) = str2num(get(H.space_RIR_mod_depth,'String'));
XStimParams.mod_phase(1) = str2num(get(H.space_RIR_mod_phase,'String'));
XStimParams.RIRpts2use = str2num(get(H.space_RIR_RIRpts2use, 'String'));
XStimParams.curr_ABL = str2num(get(H.space_RIR_ABL,'String'));
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

%%%%%%%%% check to be sure RIR FNs chosen
if isempty(FN.RIR)         %RIRs from file
    set(H.space_RIR_RIRfilepb,'value',1)
    setinfo_space_RIR
end
nFiles = length(FN.RIR);

% Add a piece of silence before and after (stim * RIR)
silence_ptsLEAD = (XStimParams.silence_lead * round(TDT.Fs/1000));
silence_ptsTRAIL = (XStimParams.silence_trail * round(TDT.Fs/1000));

%fully-cued Test
nptsNoi = round(XStimParams.curr_stimdur * TDT.Fs/1000);
Npts_totalplay = silence_ptsTRAIL + silence_ptsLEAD + nptsNoi + XStimParams.RIRpts2use + TDT.ephonefiltlen;
XStimParams.test_type = 'space_RIR FC';     str1 = 'fc';
disp('This is a FULLY CUED space_RIR test')

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

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

ISI = XStimParams.test_ISI;
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000)); 				%correct for ITD filtlength
ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));			%correct for ephonefilt

%%%%%%%Get the RIR indices
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
if get(H.space_RIR_recorddata,'Value')
    tempstr = ['    ' str1 ' space_RIR-test: ' num2str(size(XStimParams.locations,2)) ' locations;     atten: ' num2str(abs(ABL))];
    update_diary
end

Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_space_RIR_' str1 ' Temp_params;'])
clear Temp_params str1


locnum = 1:size(XStimParams.locations,2);

% update display
setInfo_space_RIR;

%Begin playing sounds
set(H.space_RIR_status,'String','Status: playing stimuli');
set(H.space_RIR_status,'BackgroundColor','green');
set(H.space_RIR_remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.space_RIRfig') & get(H.space_RIR_recorddata,'Value'))
    update_dataFN;
end

% make ramp
ramp = ones(nptsNoi,1);
ramp(1:XStimParams.ramppts) = 0:1/XStimParams.ramppts:1-(1/XStimParams.ramppts);
ramp(nptsNoi - XStimParams.ramppts+1:end) = 1-(1/XStimParams.ramppts):-1/XStimParams.ramppts:0;

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
trialnum = 1;
spikes_trial = [];
EL_trial = [];
AZ_trial = [];
repnum_trial = [];
RIR_trial = [];
Nspikes = [];
seed_trial = [];
% find indices into RIRs for each location used
nLocs = size(XStimParams.locations,2);
for iLoc = 1:nLocs
    locind(iLoc) = find(XStimParams.locations2(1,:)==XStimParams.locations(1,iLoc) & ...
        XStimParams.locations2(2,:)==XStimParams.locations(2,iLoc));
end

%%%%loop for reps
while (exist1('H.space_RIRfig') & (repnum <= XStimParams.numreps))
    tic
    %%%% loop for RIR FNs
    randseqRIRs = randperm(size(FN.RIR,2));
    for iRIR = 1:size(FN.RIR,2)
        
        randseqLocs = randperm(nLocs);
        %%%% loop for locs
        for iLoc = 1:nLocs
            %Check for pause by user
            if pause_check    return; end
            
            %%% make SAM BBN stimulus
            seed = round(rand(1)*2^30);
            noise1 = amnoise(1000,11000,XStimParams.mod_freq(1),XStimParams.mod_depth, ...
                TDT.Fs,XStimParams.curr_stimdur/1000,XStimParams.mod_phase,seed);
            noise1 = noise1 .* ramp;
            
            % normalize to ACPower
            if ~XStimParams.HiDynamicRange
                noise1 = noise1 / mom(noise1,2);
            end
            
            %Apply RIR filtering
            if pause_check    return; end        %Check for pause by user
                       
            RIR_left = RIR1{randseqRIRs(iRIR)}(locind(randseqLocs(iLoc)),:);
            RIR_right = RIR2{randseqRIRs(iRIR)}(locind(randseqLocs(iLoc)),:);
            len = length(RIR_left);
            if len < XStimParams.RIRpts2use
                RIR_left = [RIR_left zeros(1,XStimParams.RIRpts2use-len)];
                RIR_right = [RIR_right zeros(1,XStimParams.RIRpts2use-len)];
            elseif len > XStimParams.RIRpts2use
                RIR_left = RIR_left(1:XStimParams.RIRpts2use);
                RIR_right = RIR_right(1:XStimParams.RIRpts2use);
            end
                
            trial_left = conv(noise1, RIR_left);
            trial_right = conv(noise1, RIR_right);
            
            % remove DC offset
            trial_left = trial_left - mean(trial_left);
            trial_right = trial_right - mean(trial_right);
            
            % adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
            % presented through DSPs (*.eq) at 0,0 (broadband)
            trial_left = trial_left * TDT.hrtf_Lfactor;
            trial_right = trial_right * TDT.hrtf_Rfactor;
            
            %scale stimuli: note scale factor changes with XStimParams.HiDynamicRange
            trial_left = trial_left * scaleFactor;
            trial_right = trial_right * scaleFactor;
            
            %Add in the leading silent period
            trial_left =  [zeros(1,silence_ptsLEAD) trial_left];
            trial_right = [zeros(1,silence_ptsLEAD) trial_right];
            
            %Add in the trailing silent period
            trial_left =  [trial_left zeros(1,silence_ptsLAG)];
            trial_right = [trial_right zeros(1,silence_ptsLAG)];
            
            %pad with zeros for ephonefilters
            trial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
            trial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
            
            % load to buffers
            S232('push16',trial_left,Npts_totalplay);
            S232('qpop16',BUF.L1);
            S232('push16',trial_right,Npts_totalplay);
            S232('qpop16',BUF.R1);
  
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
            
            pause(Npts_totalplay/TDT.Fs +.2);
            
            while(S232('PD1status',1)) usec_delay(1000);  end
            S232('PD1stop',1);
            
            %Stop the m110 and get spiketimes
            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, round(XStimParams.curr_stimdur*2)); 			% Take 2*XStimParams.curr_stimdur spikes max
            ind = find(spikes ~= 0); 						% Get clock events that are spikes
            spikes = spikes(ind);
            ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
            if(exist1('H.space_RIRfig') & ~isempty(spikes))
                spikes = spikes(ind);
                spikes_trial = [spikes_trial; spikes/(1000/mii_us_per_sample)];
                EL_trial = [EL_trial; hrtfdirmat(1,locind(randseq(iLoc)))* ones(size(spikes))];
                AZ_trial = [AZ_trial; hrtfdirmat(2,locind(randseq(iLoc)))* ones(size(spikes))];
                repnum_trial = [repnum_trial; repnum * ones(size(spikes))];
                RIR_trial = [RIR_trial; randseqRIRs(iRIR)* ones(size(spikes))];
                seed_trial = [RIR_trial; seed* ones(size(spikes))];
                Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
            end
            ind = (iLoc-1)*nFiles + iFile;
            finalspikematrix(ind) = finalspikematrix(ind) + length(spikes);
            if pause_check    return; end
            
            remtrials = numtrials - trialnum;
            set(H.space_RIR_remtrials,'String',num2str(remtrials));
            trialnum = trialnum + 1;
            pause(0);
        end             % iLocs
        
        %Record Data
        if(exist1('H.space_RIRfig') & get(H.space_RIR_recorddata,'Value'))
            datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial RIR_trial]];
            record_data3(XStimParams,datamatrix);
        end
    end                 % iFiles
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.space_RIRfig') & ~exist1('H.space_RIR_finalspikeratefig'))
        H.space_RIR_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','space_RIR Test Spike Rate Plot',...
            'NumberTitle','off');
        H.space_RIR_spikeaxes = axes;
    end
    figure(H.space_RIR_finalspikeratefig); hold off
    %plotdiam1(XStimParams.locations, interimspikerate);
    for iLoc = 1:nLocs
        ind = [1:nFiles] + (iLoc-1)*nFiles;
        plot(1:nFiles,finalspikematrix(ind),'color',colors(rem(iLoc,10)+1,:))
        hold on
    end
    set(H.space_RIR_spikeaxes,'Color','black');
    xlabel('file#'); ylabel('spikerate'); title(['Rep # ' num2str(repnum)]);
    
    remreps = XStimParams.numreps - repnum;
    set(H.space_RIR_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.space_RIR_finalspikeratefig); hold off
set(H.space_RIR_finalspikeratefig,'Name','Final Plot for space_RIR Test');
for iLoc = 1:nLocs
    ind = [1:nFiles] + (iLoc-1)*nFiles;
    plot(1:nFiles,finalspikematrix(ind),'color',colors(rem(iLoc,10)+1,:))
    hold on
end
set(H.space_RIR_spikeaxes,'Color','black');
xlabel('file#'); ylabel('spikerate');
title('good stuff - finished!', 'FontSize',8);

set(H.space_RIR_status,'String','Status: results');
set(H.space_RIR_status,'BackgroundColor','blue');

set(H.exitspace_RIR,'Visible','on');
set(H.resetspace_RIR,'Visible','on');

% increment test number
if(exist1('H.space_RIRfig') & get(H.space_RIR_recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
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
while (exist1('H.space_RIRfig') & get(H.pausespace_RIR,'Value'))
    pause(0);
    if(~exist1('H.space_RIRfig')) return; end
    set(H.exitspace_RIR,'Visible','on');
    set(H.resetspace_RIR,'Visible','on');
    if(exist1('H.space_RIRfig') & get(H.resetspace_RIR,'Value') == 1)
        set(H.resetspace_RIR,'Value',0);
        set(H.pausespace_RIR,'Value',0);
        Reset_space_RIR;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_space_RIR;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end