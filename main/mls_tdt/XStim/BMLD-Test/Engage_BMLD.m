function [] = Engage_BMLD()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%*******************************************************************************
%	The BMLD Test operation (heavily altered from Two_source)
%**************************************************************************

rand('state',sum(100*clock));

stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);

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
    if XStimParams.space_flag
        disp('This is a FULLY-CUED BMLD test')
    else
        disp('This is an ITD BMLD test')
    end
end

%Put parameters into XStimParams
XStimParams.curr_ITD = str2num(get(H.ITD(1),'String'));
XStimParams.ITD2 = str2num(get(H.ITD(2),'String'));
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
BUF.L2				= 3;
BUF.R2				= 4;
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

%Add silence prior to stimulus
silence_len = (XStimParams.silence_lead * round(TDT.Fs/1000));
%Add silence after stimulus 
silence_len2 = (XStimParams.silence_trail * round(TDT.Fs/1000));

%Make Stimulus buffers
DUR = XStimParams.curr_stimdur;
if  XStimParams.space_flag    %fully-cued Test
    S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
else
    S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2);
    S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2);
    S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2);
    S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2);
end   

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
if XStimParams.space_flag 		% fully-cued
    S232('PD1npts',1,(silence_len2 + silence_len + DUR*(round(TDT.Fs/1000))) + TDT.ephonefiltlen + TDT.hrtffiltlen);
end

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

ITD = XStimParams.curr_ITD;
if(abs(ITD) > 250) return; end

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4atten',1,0);					% attenuation is internal to stim files
S232('PA4atten',2,0);

ISI = XStimParams.test_ISI;

%%%%%%%%%% load all HRTFs with HRTFfiletype == 2
dir = 0;
if FN.HRTFfiletype(1,2) == 2
    eval(['load -mat ' FN.space_path FN.space_std]);
    TF1_space = TF1; TF2_space = TF2;
    dir_space = dir;
end
clear dir TF1 TF2

if XStimParams.space_flag
    if FN.HRTFfiletype(1,2) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space_std]));
    elseif FN.HRTFfiletype(1,2) == 2
        hrtfdirmat = dir_space;
    end
    % find locations
    loc1 = max(find(hrtfdirmat(1,:) == XStimParams.locations(1) &...
        hrtfdirmat(2,:) == XStimParams.locations(2)));
    loc2 = max(find(hrtfdirmat(1,:) == XStimParams.locations1(1) &...
        hrtfdirmat(2,:) == XStimParams.locations1(2)));
else
    loc1 = XStimParams.curr_ITD;
    loc2 = XStimParams.ITD2;
end

% create arrays of locations and abls for all trials
nABLs = size(XStimParams.abls,2);
nLocs = XStimParams.each_loc+1;
abls = XStimParams.abls;
locations = repmat(loc2,1,nABLs);
if nLocs >1
    locations = [locations repmat(loc1,1,nABLs)];
    abls = [abls abls];
end
nTrials = size(abls,2);

%%%%%%%%%%%%%%%%%%%% make the stimuli we'll use
remreps = 1;
set(H.status,'String','Status: Building Stimuli');
set(H.remreps,'String',num2str(remreps));
irep = 1;
finalspikematrix = zeros(1,nTrials);

%Randomize the stimuli
randseq = randperm(nTrials);
itrial = 1;
while (exist1('H.BMLDfig') & (itrial <= nTrials))
    set(H.status,'BackgroundColor','red');
    %Check for pause by user
    if pause_check  return; end
    
    %Make the first stimulus
    if get(H.stim_type,'Value') ~= 9
        source1_L = get_stim(XStimParams.freq(1));
    else
        source1_L = get_stim(FN.stim_path,FN.stim);
    end
    
      % remove any DCoffset
    source1_L = source1_L - mom(source1_L,1);

    % attenuation is internal to stim files
    % scale to AC power == 1
    %source1_L = source1_L/mom(source1_L,2);

    % modulate stim1
    if ~strcmp(XStimParams.mod_type,'None')
        Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
        source1_L = source1_L .* Envelope(:)';
    end
    
    if get(H.stim_type,'Value') ~= 9      %Ramp the stimuli
        ramp_time = 5; %ms
        [source1_L] = ramp_sound(source1_L,TDT.Fs,ramp_time);
    end
    
    % remove any DCoffset
    source1_L = source1_L - mom(source1_L,1);
    
    % make R == L
    source1_R = source1_L;

    %Apply ITD filtering if conducting ITD BMLD Test
    if ~XStimParams.space_flag
        itdleft = 0; itdright = 0;
        ITD = round(loc1);
        if(ITD < 0)
            itdleft = 0;
            itdright = abs(ITD);
        elseif(ITD > 0)
            itdleft = abs(ITD);
            itdright = 0;
        end
        if(itrial == 1)
            eval(['load ' FN.ITD_path 'itdfilt' num2str(itdleft)]);
            eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
            eval(['load ' FN.ITD_path 'itdfilt' num2str(itdright)]);
            eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
        end
        source1_L = conv(source1_L,itd_filt_left);
        source1_R = conv(source1_R,itd_filt_right);
    end
    
    %Add in the leading silent period
    source1_L =  [zeros(1,silence_len) source1_L];
    source1_R = [zeros(1,silence_len) source1_R];
    
    %Add in the trailing silent period
    source1_L =  [source1_L zeros(1,silence_len2)];
    source1_R = [source1_R zeros(1,silence_len2)];
    
    %Apply HRTF filtering if conducting a virtual BMLD test
    if(XStimParams.space_flag == 1)
        if FN.HRTFfiletype(1,2) == 1
            hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locations(randseq(itrial))-1));
            hrtf_right = mtlrch([FN.space_path FN.space_std],(2*locations(randseq(itrial))));
        else
            hrtf_left = TF1_space(locations(randseq(itrial)),:);
            hrtf_right = TF2_space(locations(randseq(itrial)),:);
        end
        source1_L = conv(source1_L,hrtf_left);
        source1_R = conv(source1_R,hrtf_right);
    end
    % remove DC
    source1_L = source1_L - round(mean(source1_L));
    source1_R = source1_R - round(mean(source1_R));
    
     % these adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
    % presented through DSPs (*.eq) at 0,0 (broadband)
    source1_L = source1_L * TDT.hrtf_Lfactor;
    source1_R = source1_R * TDT.hrtf_Rfactor;    
    
    %%%%%%%%% construct the 2nd stimulus
    if get(H.stim_type2,'Value') ~= 9
        source2_L = get_stim(XStimParams.freq(2));
    else
        source2_L = get_stim(FN.stim_path2,FN.stim2);
    end
     % remove any DCoffset
     source2_L = source2_L - mom(source2_L,1);
     % attenuation is internal to stim files
     % scale to ACpower == 1
     %source2_L = source2_L/mom(source2_L,2);
    
    % modulate stim2
    if ~strcmp(XStimParams.mod_type2,'None')
        Envelope = make_env(DUR, XStimParams.mod_type2, XStimParams.mod_depth(2), XStimParams.mod_freq(2), XStimParams.mod_phase(2));
        source2_L = source2_L .* Envelope(:)';
    end
    
    if get(H.stim_type2,'Value') ~= 9  	%Ramp the stimuli
        ramp_time = 5; %ms
        source2_L = ramp_sound(source2_L,TDT.Fs,ramp_time);
    end
    
    % remove any DCoffset
    source2_L = source2_L - mom(source2_L,1);
    
        % start with R == L
    source2_R = source2_L;

    %Apply ITD filtering if conducting ITD-based BMLD Test
    if ~XStimParams.space_flag
        itdleft = 0; itdright = 0;
        ITD = round(locations(randseq(itrial)));
        if(ITD < 0)
            itdleft = 0;
            itdright = abs(ITD);
        elseif(ITD > 0)
            itdleft = abs(ITD);
            itdright = 0;
        end
        if(itrial == 1)
            eval(['load ' FN.ITD_path 'itdfilt' num2str(itdleft)]);
            eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
            eval(['load ' FN.ITD_path 'itdfilt' num2str(itdright)]);
            eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
        end
        source2_L = conv(source2_L,itd_filt_left);
        source2_R = conv(source2_R,itd_filt_right);
    end
    
    %Add in the leading silent period
    source2_L =  [zeros(1,silence_len) source2_L];
    source2_R = [zeros(1,silence_len) source2_R];
    
    %Add in the trailing silent period
    source2_L =  [source2_L zeros(1,silence_len2)];
    source2_R = [source2_R zeros(1,silence_len2)];
    
    %Apply HRTF filtering if virtual BMLD test
    if XStimParams.space_flag
        if FN.HRTFfiletype(1,2) == 1
            hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locations(randseq(itrial))-1));
            hrtf_right = mtlrch([FN.space_path FN.space_std],2*locations(randseq(itrial)));
        else
            hrtf_left = TF1_space(locations(randseq(itrial)),:);
            hrtf_right = TF2_space(locations(randseq(itrial)),:);
        end
        source2_L = conv(source2_L,hrtf_left);
        source2_R = conv(source2_R,hrtf_right);
    end
    
    % remove DC offset
    source2_L = source2_L - round(mean(source2_L));
    source2_R = source2_R - round(mean(source2_R));

    % these adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
    % presented through DSPs (*.eq) at 0,0 (broadband)
    source2_L = source2_L * TDT.hrtf_Lfactor;
    source2_R = source2_R * TDT.hrtf_Rfactor;    
    
    trial_left = (source1_L * XStimParams.factor(1) + source2_L * XStimParams.factor(2))/2;
    trial_right = (source1_R * XStimParams.factor(1) + source2_R * XStimParams.factor(2))/2;
    
    % 3/7/07
    trial_left = trial_left * TDT.scaleFactor;
    trial_right = trial_right * TDT.scaleFactor;
    
    %pad with zeros
    filttrial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
    filttrial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
    
    % save stims to disk with name of locations _ abls
    if(exist1('H.BMLDfig'));
        S232('push16',filttrial_left,length(filttrial_left));
        S232('qpop16',BUF.L1);
        fname = ['stimbuf_left_' num2str(locations(randseq(itrial))) ...
                '_' num2str(abls(randseq(itrial)))];
        evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
                [' ''' stimuli_dir fname ''' ']   ',0);'];
        eval(evalstr);
        temp_left = dama2pc(BUF.L1);
        S232('push16',filttrial_right,length(filttrial_right));
        S232('qpop16',BUF.R1);
        fname = ['stimbuf_right_' num2str(locations(randseq(itrial))) ...
                '_' num2str(abls(randseq(itrial)))];
        evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
                [' ''' stimuli_dir fname ''' ']   ',0);'];
        eval(evalstr);
        temp_right = dama2pc(BUF.R1);
    end
    
    remtrials = nTrials - itrial;
    set(H.remtrials,'String',num2str(remtrials));
    itrial = itrial + 1;
    set(H.status,'BackgroundColor','blue');
    pause(0);
end 										%end loop over trials

%%%%%%%%%%%%%%%%%%%%%%%%%%% finished making sounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin playing sounds   
set(H.status,'String','Status: Playing Stimuli');
set(H.status,'BackgroundColor','green');
set(H.remreps,'String',num2str(XStimParams.numreps));
irep = 1;
datamatrix = [];
rep_spikes = [];

% increment testnumber
if(exist1('H.BMLDfig') & get(H.recorddata,'Value'))
    update_dataFN;
end

%loop for reps
while (exist1('H.BMLDfig') & (irep <= XStimParams.numreps))
    %Randomize the stimuli
    randseq = randperm(nTrials);   
    itrial = 1;
    spikes_trial = [];
    ABL_trial = [];
    EL_trial = [];
    AZ_trial = [];
    irep_trial = [];
    Nspikes = [];
    rep_spikes = [rep_spikes; zeros(1,nTrials)];                  % n spikes for this rep
    
    % loop for trials
    tic
    while (exist1('H.BMLDfig') & (itrial <= nTrials+1))
        
        %Check for pause by user
        if pause_check  return; end
        
        %Wait till PD1 is finished
        while S232('PD1status',1) usec_delay(1000); end
        
        if(itrial <= nTrials)
            fname = ['stimbuf_left_' num2str(locations(randseq(itrial))) ...
                    '_' num2str(abls(randseq(itrial)))];
            evalstr = ['S232(''disk2dama16'',BUF.L1,'  [' ''' stimuli_dir fname ''' '] ',0);'];
            eval(evalstr);
            fname = ['stimbuf_right_' num2str(locations(randseq(itrial))) ...
                    '_' num2str(abls(randseq(itrial)))];
            evalstr = ['S232(''disk2dama16'',BUF.R1,' [' ''' stimuli_dir fname ''' '] ',0);'];
            eval(evalstr);
        end
        
        %Wait till PD1 is finished
        while S232('PD1status',1) usec_delay(1000); end
        S232('PD1stop',1);
        
        %Stop the m110 and get spikes
        if(itrial > 1)							% first trial just for loading
            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, 1000); 			% Take 100 spikes max
            ind = find(spikes ~= 0); 						% Get clock events that are spikes
            spikes = spikes(ind);
            ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
            if(exist1('H.BMLDfig') & ~isempty(spikes)) 
                spikes = spikes(ind);
                spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
                EL_trial = [EL_trial;locations(randseq(itrial-1))* ones(size(spikes))];
                AZ_trial = [AZ_trial;locations(randseq(itrial-1))* ones(size(spikes))];
                ABL_trial = [ABL_trial;abls(randseq(itrial-1))* ones(size(spikes))];
                irep_trial = [irep_trial;irep * ones(size(spikes))];
                Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
                rep_spikes(irep,randseq(itrial-1)) = length(spikes);
            end
        end
        
        %Check for pause by user
        if pause_check  return; end
        
        if(exist1('H.BMLDfig') & (itrial <= nTrials))
            S232('seqplay',BUF.playspec1);
            S232('PD1arm',1);
            %Send trigger
            %Set up MII
            m100x( C_.INIT );
            m110dx( C_.INIT );
            m110dx( C_.CLOCK, mii_us_per_sample);
            m110dx( C_.MODE, M_.PST );
            if (itrial <= nTrials)
                while toc < ISI/1000     end
                %Start clock
                m110dx( C_.START);
                %Send pulse: PD1 GO!
                m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
                tic
            end
        end
        
        remtrials = nTrials - itrial +1;
        set(H.remtrials,'String',num2str(remtrials));
        itrial = itrial + 1;
        pause(0);
    end %end loop over trials
    
    %Plot BMLD plot
    if exist1('H.BMLDfig') & ~exist1('H.BMLDspikefig')
        H.BMLDspikefig = figure('Position',[700 20 550 500],...
            'Name','BMLD Test Spike Rate Plot',...
            'NumberTitle','off');
        H.spikeaxes = axes;
        hold on
        set(H.spikeaxes,'Color','black');
        xlabel('relative attenuation (dB)'); ylabel('spikes');
    end
    figure(H.BMLDspikefig)
    title(['after ' num2str(irep) '  reps']);
    plot(abls(1:nABLs),rep_spikes(irep,1:nABLs),'m.')
    if nLocs == 2
        plot(abls(nABLs+1:end),rep_spikes(irep,nABLs+1:end),'y.')
    end        
    %Record Data
    if(exist1('H.BMLDfig') & get(H.recorddata,'Value'))
        tempseq{irep} = randseq;
        datamatrix = [datamatrix;[Nspikes spikes_trial irep_trial EL_trial AZ_trial ABL_trial]];
        record_data3(XStimParams,datamatrix, tempseq);
    end
    
    remreps = XStimParams.numreps - irep;
    set(H.remreps,'String',num2str(remreps));
    irep = irep + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure   BMLD ????
figure(H.BMLDspikefig)
title('BMLD final')
plot(abls(1:nABLs),mean(rep_spikes(:,1:nABLs),1),'r*')
plot([abls(1:nABLs); abls(1:nABLs)],[mean(rep_spikes(:,1:nABLs))+ ...
    std(rep_spikes(:,1:nABLs)); mean(rep_spikes(:,1:nABLs))-std(rep_spikes(:,1:nABLs))],'r')

if nLocs == 2
    plot(abls(nABLs+1:end),mean(rep_spikes(:,nABLs+1:end),1),'g*')
    plot([abls(nABLs+1:end); abls(1:nABLs)],[mean(rep_spikes(:,nABLs+1:end))+ ...
        std(rep_spikes(:,nABLs+1:end));  mean(rep_spikes(:,nABLs+1:end))-std(rep_spikes(:,nABLs+1:end))],'g')
end        

ind = find(abls == min1(abls));
noise_alone = rep_spikes(:,ind); noise_alone = noise_alone(:);
plot(abls(1:nABLs),ones(1,nABLs)*mean(noise_alone),'b*')
plot([min(abls) max(abls)],[mean(noise_alone)+ std(noise_alone) ...
        mean(noise_alone)+ std(noise_alone)],'b-')
plot([min(abls) max(abls)],[mean(noise_alone)- std(noise_alone) ...
        mean(noise_alone)- std(noise_alone)],'b-')
clear ind


set(H.status,'BackgroundColor','blue');
set(H.status,'String','Status: Results');
set(H.exitBMLD,'Visible','on');
set(H.resetBMLD,'Visible','on');

% increment test number
if(exist1('H.BMLDfig') & get(H.recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end
S232('PA4mute',1);
S232('PA4mute',2);

%%%%%%%%%
function [stim] = get_stim(param1, param2)
% param1: XStimParams.freq or FN.stim_path
% param2: FN.stim

global H
global XStimParams
global TDT

switch get(H.stim_type,'Value')      
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
        disp('Stimulus type not supported for BMLD tests.  Reset to BROADBAND');
        return
end
stim = stim / mom(stim,2);

%%%%%%
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
while (exist1('H.BMLDfig') & get(H.pauseBMLD,'Value'))
    pause(0);
    if(~exist1('H.BMLDfig')) return; end         
    set(H.exitBMLD,'Visible','on');
    set(H.resetBMLD,'Visible','on');
    if(exist1('H.BMLDfig') & get(H.resetBMLD,'Value') == 1)
        set(H.resetBMLD,'Value',0);
        set(H.pauseBMLD,'Value',0);
        Reset_BMLD;   flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_BMLD;   flag=1;
        return;
    end
end

if XStimParams.reset_flag
    flag = 1;
    XStimParams.reset_flag = 0;
end