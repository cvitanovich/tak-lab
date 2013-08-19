function [] = Engage_Delay()

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


% //////////
% IMPORTANT: Revert to random sequence!!!!!!!!!!!!!!! 
UseRandPlaySeq = 1;

SaveStimToDisk = 0;  % 1 = save one repetiion to disk and play, 0 = Generate stimuli on the fly

plotAllStim = 0; % plot all stim as they are generated
plotStimTrial = 0; % plot this trial's stim.

% XStimParams added to XStimParams_func  Run to add!!!
%    'SilentDelay',              0, ...
%    'uncorrel',                 0, ...
%    'SilentLag',                0, ...
%    'randOnsetPerms',           0, ...
%    'DelayMods',                [0 0 0 0], ...
%    'DelayTimes',               [-200 -100 -50 -20 -10 -5 -2 -1 1 2 5 10 20 50 100 200], ...
%    'ramp_timeS_on',            [1 2.5 5 10 20 50], ...
%    'ramp_timeS_off',           [1 2.5 5 10 20 50], ...
%    'lag_seg_ext',              [0 30 NaN], ...
%    'OnOff_mode',               1, ...
%    'Delay_default',            1 ...

%     XStimParams.SilentDelay =0;
%     XStimParams.uncorrel =0;
%     XStimParams.SilentLag =0;
%     XStimParams.randOnsetPerms =0;
%     XStimParams.DelayMods = [1 2 3 4];
%     XStimParams.DelayTimes = [-200 -100 -50 -20 -10 -5 -2 -1 1 2 5 10 20 50 100 200];
%     XStimParams.ramp_timeS_on = [1 2.5 5 10 20 50];
%     XStimParams.ramp_timeS_off = [1 2.5 5 10 20 50];
%     XStimParams.lag_seg_ext = [0 30 NaN];
%     XStimParams.OnOff_mode = 1;
%     XStimParams.Delay_default = 1;

%if get(H.Delay_stim_type,'Value')==2; %else == 2
    if PEstim(0)== 0 
        warning('Random numbers have changed. Use stimulus files');
        return;
    end
%end

% get stimulus parameters from Delay figure
if(XStimParams.randOnsetPerms)
    % Get rid of NaN delays
    DelayTimes=zeros(1,1);
    mp=1;
    for mt = 1:size(XStimParams.DelayTimes,2)
        if(~isnan(XStimParams.DelayTimes(mt)))
            DelayTimes(mp)=XStimParams.DelayTimes(mt);
            mp=mp+1;
        end
    end

    MaxDelay_ms=max1(abs(DelayTimes));
    MaxDelay_pnts=round(TDT.Fs*(MaxDelay_ms/1000));

    
    % Get rid of NaN onset ramps
    ramp_timeS_on=zeros(1,1);
    mp=1;
    for mt = 1:size(XStimParams.ramp_timeS_on,2)
        if(~isnan(XStimParams.ramp_timeS_on(mt)))
            ramp_timeS_on(mp)=XStimParams.ramp_timeS_on(mt);
            mp=mp+1;
        end
    end
    MaxOnset=max1(abs(ramp_timeS_on));
    
    % Get rid of NaN offset ramps
    ramp_timeS_off=zeros(1,1);
    mp=1;
    for mt = 1:size(XStimParams.ramp_timeS_off,2)
        if(~isnan(XStimParams.ramp_timeS_off(mt)))
            ramp_timeS_off(mp)=XStimParams.ramp_timeS_off(mt);
            mp=mp+1;
        end
    end
    MaxOffset=max1(abs(ramp_timeS_off)); % ramp_timeS_off ???????
    
    lag_seg_ext=zeros(1,1);
    mp=1;
    for mt = 1:size(XStimParams.lag_seg_ext,2)
        if(~isnan(XStimParams.lag_seg_ext(mt)))
            lag_seg_ext(mp)=XStimParams.lag_seg_ext(mt);
            mp=mp+1;
        end
    end
    MaxLag_ms=(max1(abs(lag_seg_ext)));  
      
else
    SilentOnset=XStimParams.SilentDelay;
    if(SilentOnset & XStimParams.DelayTimes(1) >= XStimParams.curr_stimdur)
        XStimParams.DelayTimes(1)=XStimParams.curr_stimdur-10;
        set(H.delaytime,'String', num2str(XStimParams.DelayTimes(1)));
        warndlg('delay set to stim. dur. minus 10 ms', 'warning!');
        waitforbuttonpress;
    end
    Delay_pnts=round(TDT.Fs*(XStimParams.DelayTimes(1)/1000));
    uncorrel=XStimParams.uncorrel;
    SilentLag=XStimParams.SilentLag;
    MaxDelay_ms=abs(XStimParams.DelayTimes(1));
    MaxDelay_pnts=abs(Delay_pnts);
    
    ramp_time_on = XStimParams.ramp_timeS_on(1); %ms
    ramp_time_off = XStimParams.ramp_timeS_off(1); %ms
    MaxLag_ms = abs(XStimParams.lag_seg_ext(1));
end
%MaxLag_Pnts = abs(MaxLag_ms*round(TDT.Fs/1000));
MaxLag_Pnts=round(TDT.Fs*(MaxLag_ms/1000));

%disp([MaxDelay_pnts MaxLag_Pnts]);


% 50 random noise states (corresponding with saved files).
RandStates = [...
73 167 617 364 662 593 538 194 853 610 294 ...
479 71 105 162 770 143 116 252 101 377 706 ...
273 574 915 661 935 355 530 540 220 232 886 ...
70 65 571 35 339 87 281 795 283 974 248 ...
995 936 769 943 127 224];

if(findobj('Tag','Raster'))
    close('Raster');
end

% New state for generating trials
SeqRandState=sum(100*clock);
rand('state',SeqRandState);
 
stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);

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

%Put parameters into XStimParams
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.reset_flag = 0;

% if (get(H.stim_type,'Value') == 9) & isempty(FN.stim)
%     [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select stimulus File');
% end

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

%Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
silence_len = (XStimParams.silence_lead * round(TDT.Fs/1000));
%Add a piece of silence after stimulus 
silence_len2 = (XStimParams.silence_trail * round(TDT.Fs/1000));

%Make Stimulus buffers
DUR = XStimParams.curr_stimdur;

%if XStimParams.space_flag == 1	    %fully-cued Test
    S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen + abs(MaxDelay_pnts) + MaxLag_Pnts);
    S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen + abs(MaxDelay_pnts) + MaxLag_Pnts);
    S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen + abs(MaxDelay_pnts) + MaxLag_Pnts);
    S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen + abs(MaxDelay_pnts) + MaxLag_Pnts);
    %end   

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
%if XStimParams.space_flag == 1	% fully-cued
    S232('PD1npts',1,(silence_len2 + silence_len + DUR*(round(TDT.Fs/1000))) + TDT.ephonefiltlen + TDT.hrtffiltlen + abs(MaxDelay_pnts) + MaxLag_Pnts);
%end


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
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

ISI = XStimParams.test_ISI;

%Get all HRTF spectra indices
XStimParams.locations = GUI.locations1';
if ~size(XStimParams.locations,2) % | size(XStimParams.locations,2)>1 % Limit trials to one location
    set(H.pickerfig,'Color', [.1 .2 .8]);
    set(H.picker_error,'visible','on');
    pause;
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end


%%%%%%%%%% load all HRTFs with HRTFfiletype == 2
  dir = 0;
if FN.HRTFfiletype(1) == 2
    eval(['load -mat ' FN.space_path FN.space]);
    TF1_space = TF1; TF2_space = TF2;
    dir_space = dir;
end
if FN.HRTFfiletype(2) == 2
    eval(['load -mat ' FN.ILA_path FN.ildalone]);
    TF1_ila = TF1; TF2_ila = TF2;
    dir_ILA = dir;
end
if FN.HRTFfiletype(3) == 2
    eval(['load -mat ' FN.ITA_path FN.itdalone]);
    TF1_ita = TF1; TF2_ita = TF2;
    dir_ITA = dir;
end
clear dir TF1 TF2

%if XStimParams.space_flag == 1
    if FN.HRTFfiletype(1) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space]));
    elseif FN.HRTFfiletype(1) == 2
        hrtfdirmat = dir_space;
    end
%end
%%%%%%%%%%%

% find locations for 2 or 3 sources
clear locind locations
numlocs = 0;

offset_el = XStimParams.offset_el(2);%/2;
offset_az = XStimParams.offset_az(2);%/2;

for locnum = 1:size(XStimParams.locations,2)
    loc1 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,locnum)) &...
        hrtfdirmat(2,:) == (XStimParams.locations(2,locnum))));
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

%get reference Lref, Rref and ABL for ABLalone test (always use just first ref location)
if XStimParams.ABLalone_flag 
    ind00 = max(find(hrtfdirmat(1,:) == XStimParams.el & hrtfdirmat(2,:) == XStimParams.az));
    if FN.HRTFfiletype(1) == 1
        Lref = mtlrch([FN.space_path FN.space],(2*ind00)-1);
        Rref = mtlrch([FN.space_path FN.space],2*ind00);
    else
        Lref = TF1_space(ind00,:)';
        Rref = TF2_space(ind00,:)';
    end
    tempL = ERBFilterBankB(Lref, fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
    tempR = ERBFilterBankB(Rref, fcoefs) .* Factormat;
    [ildref ablref] = calclevel_time(tempL,tempR, cF);
    %ITDref = calcitd(tempL,tempR, cF, Fs, ones(n_cF,1));
    clear tempL tempR temp
end

%%%%%%%%%%%%%%%%%%%% make the stimuli we'll use
remreps = 1;
set(H.status,'String','Status: Building Stimuli');
set(H.remreps,'String',num2str(remreps));
repnum = 1;

% generate stimulus sequences
if(XStimParams.randOnsetPerms==0) 
    % sequence spans locations only
    numtrials = size(XStimParams.locations,2);
    LocationsSeq=locations;
    LocIndSeq=locind;
    
    DelayModsSeq= zeros(numtrials,1);
    %  = [1 2 3]; % 1=normal, 2=uncorr, 3=corr without onset (lead segment), 4 corr without offset (lag segment)
    for tr = 1:numtrials
        if(uncorrel & ~SilentOnset & ~SilentLag)
            DelayModsSeq(tr)=2;
        else
            if(~uncorrel & SilentOnset & ~SilentLag)
                DelayModsSeq(tr)=3;
            else
                if(~uncorrel & ~SilentOnset & SilentLag)
                    DelayModsSeq(tr)=4;
                else
                    DelayModsSeq(tr)=1;
                end
            end
        end
    end
    
        clear DelayTimesSeq
    DelayTimesSeq = zeros(numtrials,1)+XStimParams.DelayTimes(1);
    %disp(DelayTimesSeq);
    
        clear OnsetSeq
    OnsetSeq = zeros(numtrials,1)+XStimParams.ramp_timeS_on(1);
        clear OffsetSeq
    OffsetSeq = zeros(numtrials,1)+XStimParams.ramp_timeS_off(1);
    
        clear LagExtSegsSeq
    LagExtSegsSeq = zeros(numtrials,1)+XStimParams.lag_seg_ext(1);
    
        clear RandStateNoiseSeq
    StateLow = 1; StateHi = 50; % choose state between 1 and 50
    RandStateNoiseSeq = zeros(numtrials,1)+RandStates(round(StateLow + (StateHi-StateLow) * rand(1)));
    
    
else % random parameters
    % getting these here because they don't seem to be updated properly...
    XStimParams.uncorrel = get(H.uncorrel,'Value');
    uncorrel=XStimParams.uncorrel;
    XStimParams.SilentDelay = get(H.silentdelay,'Value');
    SilentOnset=XStimParams.SilentDelay;
    XStimParams.SilentLag = get(H.silentlag,'Value');
    SilentLag=XStimParams.SilentLag;
    
    if(uncorrel & SilentOnset & SilentLag)
        XStimParams.DelayMods = [1 2 3 4]; 
        % 1=normal, 2=uncorr, 3=corr no lead segment, 4=corr no lag segment
    else
        if(uncorrel & SilentOnset & ~SilentLag)
            %XStimParams.DelayMods = [1 2 3];
            XStimParams.DelayMods = [1 2 3 5]; 
            % 1=normal, 2=uncorr, 3=corr no lead segment, , 5=uncorr & no lead segment
        else
            if(~uncorrel & SilentOnset & ~SilentLag)
                XStimParams.DelayMods = [1 3]; % 3=corr without onset
            else
                if(uncorrel & ~SilentOnset & ~SilentLag)
                    XStimParams.DelayMods = [1 2]; % 1=normal, 2=uncorr
                else
                    if(~uncorrel & ~SilentOnset & SilentLag)
                        XStimParams.DelayMods = [1 4]; % 1=normal, 4=corr no lag segment
                    else
                       if(~uncorrel & ~SilentOnset & ~SilentLag)
                            XStimParams.DelayMods = [1]; % 1=normal
                       else
                           warning('no valid stimulus combination selected');
                           XStimParams.DelayMods = [1]; % 1=normal
                       end
                    end
                end
            end
        end
    end
     
    % sequence spans  locations, delays, and onset mods, on/off ramps
    numDelayMods = size(XStimParams.DelayMods,2);
    numDelays = size(DelayTimes,2);
    numlocations = size(locations,2);
    numOnRamps = size(ramp_timeS_on,2);
    numOffRamps = size(ramp_timeS_off,2);
    numLagExtSegs = size(lag_seg_ext,2);
    
    numtrials=numlocations*numDelays*numOnRamps*numOffRamps*numLagExtSegs*numDelayMods;
    
    clear LocationsSeq LocIndSeq
    % set locations in stim sequence
    LocationsSeq= zeros(2,numtrials); % az & el
    LocIndSeq= zeros(2,numtrials);
    perm=1;
    for tr = 1:numtrials
        if(tr >= (numtrials/numlocations)*perm)
            LocationsSeq(:,tr)=locations(:,perm);
            LocIndSeq(:,tr)=locind(:,perm);
            perm=perm+1;
        else
            LocationsSeq(:,tr)=locations(:,perm);
            LocIndSeq(:,tr)=locind(:,perm);
        end
    end
    
    clear DelayTimesSeq
    % set lag delays in stim sequence
    DelayTimesSeq= zeros(numtrials,1);
    perm=1;
    subp=1;
    for tr = 1:numtrials
        if(tr >= ((numtrials/numlocations)/numDelays)*perm)
            DelayTimesSeq(tr)=DelayTimes(subp);
            perm=perm+1;
            subp=subp+1;
            if(subp>numDelays)
                subp=1;
            end
        else
            DelayTimesSeq(tr)=DelayTimes(subp);
        end
    end
    
        clear OnsetSeq
    % set onset mods in stim sequence
    OnsetSeq= zeros(numtrials,1);
    perm=1;
    subp=1;
    for tr = 1:numtrials
        if(tr >= (((numtrials/numlocations)/numDelays)/numOnRamps)*perm)
            OnsetSeq(tr)=ramp_timeS_on(subp);
            perm=perm+1;
            subp=subp+1;
            if(subp>numOnRamps)
                subp=1;
            end
        else
            OnsetSeq(tr)=ramp_timeS_on(subp);
        end
    end
    
        clear OffsetSeq
    % set onset mods in stim sequence
    OffsetSeq= zeros(numtrials,1);
    perm=1;
    subp=1;
    for tr = 1:numtrials
        if(tr >= ((((numtrials/numlocations)/numDelays)/numOnRamps)/numOffRamps)*perm)
            OffsetSeq(tr)=ramp_timeS_off(subp);
            perm=perm+1;
            subp=subp+1;
            if(subp>numOffRamps)
                subp=1;
            end
        else
            OffsetSeq(tr)=ramp_timeS_off(subp);
        end
    end
    
    clear LagExtSegsSeq
    % set Lag Extension Segment mods in stim sequence
    LagExtSegsSeq= zeros(numtrials,1);
    perm=1;
    subp=1;
    for tr = 1:numtrials
        if(tr >= (((((numtrials/numlocations)/numDelays)/numOnRamps)/numOffRamps)/numLagExtSegs )*perm)
            LagExtSegsSeq(tr)=lag_seg_ext(subp);
            perm=perm+1;
            subp=subp+1;
            if(subp>numLagExtSegs)
                subp=1;
            end
        else
            LagExtSegsSeq(tr)=lag_seg_ext(subp);
        end
    end
    
        clear DelayModsSeq
    % set onset mods in stim sequence
    DelayModsSeq= zeros(numtrials,1);
    perm=1;
    subp=1;
    for tr = 1:numtrials
        if(tr >= ((((((numtrials/numlocations)/numDelays)/numOnRamps)/numOffRamps)/numLagExtSegs)/numDelayMods )*perm)
            DelayModsSeq(tr)=XStimParams.DelayMods(subp);
            perm=perm+1;
            subp=subp+1;
            if(subp>numDelayMods)
                subp=1;
            end
        else
            DelayModsSeq(tr)=XStimParams.DelayMods(subp);
        end
    end
    
    % make INITIAL noise state sequence  RandStates
    RandStateNoiseSeq= zeros(numtrials,1);
    tr=1;
    StateLow = 1; StateHi = 50; % choose state between 1 and 50
    while tr <= size(RandStateNoiseSeq,1);
        RandStateNoiseSeq(tr)=RandStates(round(StateLow + (StateHi-StateLow) * rand(1)));
        tr=tr+1;
    end
    
    % ///////////////////////////
    % Eliminate unnecessary trials
    
    % Eliminate trials where delay exceeds stim. duration - on/offset
    % and where the onset or offset segment is silenced
    OrigNumtrials=numtrials;
    numtrials=1;
    tr=1;
    while tr <= size(DelayModsSeq,1);
        if((DelayModsSeq(tr)==3 | DelayModsSeq(tr)==4 | DelayModsSeq(tr)==5) & abs(DelayTimesSeq(tr)) > XStimParams.curr_stimdur-(MaxOnset+MaxOffset) )
            LocationsSeq(:,tr)=[,];
            LocIndSeq(:,tr)=[,];
            DelayTimesSeq(tr)=[];
            OnsetSeq(tr)=[];
            OffsetSeq(tr)=[];
            LagExtSegsSeq(tr)=[];
            DelayModsSeq(tr)=[];
            RandStateNoiseSeq(tr)=[];
        else
            tr=tr+1;
            numtrials=numtrials+1;
        end
    end
    numtrials=numtrials-1;
    if(numtrials<OrigNumtrials)
        wrntext = ['Eliminated ' num2str(OrigNumtrials-numtrials) ' trials where delay exceeded stim. dur. - on/offset (' num2str((MaxOnset+MaxOffset)) ' ms) but was silenced'];
        disp(wrntext);
    end
    
    % Eliminate trials where the lag segment is extended (DelayModsSeq == 4) but is silenced
    OrigNumtrials=numtrials;
    numtrials=1;
    tr=1;
    while tr <= size(DelayModsSeq,1);
        if(DelayModsSeq(tr)==4 & LagExtSegsSeq(tr) > 0 )
            LocationsSeq(:,tr)=[,];
            LocIndSeq(:,tr)=[,];
            DelayTimesSeq(tr)=[];
            OnsetSeq(tr)=[];
            OffsetSeq(tr)=[];
            LagExtSegsSeq(tr)=[];
            DelayModsSeq(tr)=[];
            RandStateNoiseSeq(tr)=[];
        else
            tr=tr+1;
            numtrials=numtrials+1;
        end
    end
    numtrials=numtrials-1;
    if(numtrials<OrigNumtrials)
        wrntext = ['Eliminated ' num2str(OrigNumtrials-numtrials) ' trials where the lag segment is both extended and silenced'];
        disp(wrntext);
    end
    
    
    % Eliminate trials where DelayTime < -DUR but the lag is extended or silenced
%     % This might now be accomplished above???????????
%     OrigNumtrials=numtrials;
%     numtrials=1;
%     tr=1;
%     while tr <= size(DelayTimesSeq,1);
%         if(DelayTimesSeq(tr) < -DUR & (LagExtSegsSeq(tr) > 0 | DelayModsSeq(tr) == 4) )
%             LocationsSeq(:,tr)=[,];
%             LocIndSeq(:,tr)=[,];
%             DelayTimesSeq(tr)=[];
%             OnsetSeq(tr)=[];
%             OffsetSeq(tr)=[];
%             LagExtSegsSeq(tr)=[];
%             DelayModsSeq(tr)=[];
%             RandStateNoiseSeq(tr)=[];
%         else
%             tr=tr+1;
%             numtrials=numtrials+1;
%         end
%     end
%     numtrials=numtrials-1;
%     if(numtrials<OrigNumtrials)
%         wrntext = ['Eliminated ' num2str(OrigNumtrials-numtrials) ' trials where a negative lag delay exceedes stim duration but is extended or silenced'];
%         disp(wrntext);
%     end
        
    % Eliminate trials where lag is extended and  delay > 10 | < 5 ms
    %                                                     ///////////
    OrigNumtrials=numtrials;
    numtrials=1;
    tr=1;
    while tr <= size(DelayTimesSeq,1);
        if( (abs(DelayTimesSeq(tr)) > 10 |abs(DelayTimesSeq(tr)) < 5) & LagExtSegsSeq(tr) > 0)
        
            LocationsSeq(:,tr)=[,];
            LocIndSeq(:,tr)=[,];
            DelayTimesSeq(tr)=[];
            OnsetSeq(tr)=[];
            OffsetSeq(tr)=[];
            LagExtSegsSeq(tr)=[];
            DelayModsSeq(tr)=[];
            RandStateNoiseSeq(tr)=[];
        else
            tr=tr+1;
            numtrials=numtrials+1;
        end
    end
    numtrials=numtrials-1;
    if(numtrials<OrigNumtrials)
        wrntext = ['Eliminated ' num2str(OrigNumtrials-numtrials) ' trials where lag is extended and  delay > 10 | < 5 ms'];
        disp(wrntext);
    end
    
    % Eliminate trials where the lag extenstion is negative but the lag segment is >= 15 ms
    %                                                                      ///////////
    OrigNumtrials=numtrials;
    numtrials=1;
    tr=1;
    while tr <= size(DelayTimesSeq,1);
        if( abs(DelayTimesSeq(tr)) + LagExtSegsSeq(tr) >= 15  & LagExtSegsSeq(tr) < 0)
        
            LocationsSeq(:,tr)=[,];
            LocIndSeq(:,tr)=[,];
            DelayTimesSeq(tr)=[];
            OnsetSeq(tr)=[];
            OffsetSeq(tr)=[];
            LagExtSegsSeq(tr)=[];
            DelayModsSeq(tr)=[];
            RandStateNoiseSeq(tr)=[];
        else
            tr=tr+1;
            numtrials=numtrials+1;
        end
    end
    numtrials=numtrials-1;
    if(numtrials<OrigNumtrials)
        wrntext = ['Eliminated ' num2str(OrigNumtrials-numtrials) ' trials where the lag extenstion is negative but the lag segment is >= 15 ms'];
        disp(wrntext);
    end 
    
    % Eliminate trials where NEGATIVE lag Extensions are longer than or equal to lag delay
    OrigNumtrials=numtrials;
    numtrials=1;
    tr=1;
    while tr <= size(DelayTimesSeq,1);
        if( abs(DelayTimesSeq(tr)) <= -LagExtSegsSeq(tr) & ...
            LagExtSegsSeq(tr) < 0 )
        
            LocationsSeq(:,tr)=[,];
            LocIndSeq(:,tr)=[,];
            DelayTimesSeq(tr)=[];
            OnsetSeq(tr)=[];
            OffsetSeq(tr)=[];
            LagExtSegsSeq(tr)=[];
            DelayModsSeq(tr)=[];
            RandStateNoiseSeq(tr)=[];
        else
            tr=tr+1;
            numtrials=numtrials+1;
        end
    end
    numtrials=numtrials-1;
    if(numtrials<OrigNumtrials)
        wrntext = ['Eliminated ' num2str(OrigNumtrials-numtrials) ' trials where NEGATIVE lag Extensions are longer than or equal to lag delay'];
        disp(wrntext);
    end
    
end % random parameters

disp(['trial ',...
'elevation ',...
'azimuth ',...
'delay ',...
'onset ',...
'offset ',...
'LagSegExt ',...
'mod (1=async 2=uncorr 3=no lead seg 4=no lag seg)'...
]);

SeqOrder = zeros(9,numtrials);
for d = numtrials:-1:1
    SeqOrder(1,d) = d;
    SeqOrder(2,d) = LocationsSeq(1,d);
    SeqOrder(3,d) = LocationsSeq(2,d);
    SeqOrder(4,d) = DelayTimesSeq(d);
    SeqOrder(5,d) = OnsetSeq(d);
    SeqOrder(6,d) = OffsetSeq(d);
    SeqOrder(7,d) = LagExtSegsSeq(d);
    SeqOrder(8,d) = DelayModsSeq(d);
    %SeqOrder(9,d) = RandStateNoiseSeq(d);
end
for d = numtrials:-1:1
    disp([num2str(SeqOrder(1,d)), ' ',...
        num2str(SeqOrder(2,d)), ' ',...
        num2str(SeqOrder(3,d)), ' ',...
        num2str(SeqOrder(4,d)), ' ',...
        num2str(SeqOrder(5,d)), ' ',...
        num2str(SeqOrder(6,d)), ' ',...
        num2str(SeqOrder(7,d)), ' ',...
        num2str(SeqOrder(8,d)), ' ',...
        %num2str(SeqOrder(9,d)), ' ',...
        ]);
end
clear('SeqOrder');


finalspikematrix = zeros(1,numtrials);
StimDur = DUR + MaxLag_ms + MaxDelay_ms;

if(SaveStimToDisk) % Save stimuli to disk?
	%repnum = 1;
	%Randomize the stimuli
	%randseq = randperm(numtrials);
	%disp(randseq);
	randseq = [1:numtrials]; % randomized when playing
	trialnum = 1;
	%loop for reps
	%while (exist1('H.Delayfig') & (repnum <= XStimParams.numreps))
		while (exist1('H.Delayfig') & (trialnum <= numtrials))
            set(H.status,'BackgroundColor','red');
            %Check for pause by user
            if pause_check  return; end
                  
            % ////////////////////////////////////////////
            %Make the stimuli
            useFile=0;
            if get(H.Delay_stim_type,'Value') == 2  % Generate
                useFile=1;
            end
            whichRamp=get(H.DelayOnOff_mode,'Value'); % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
            
            %StimDur = DUR + MaxLag_ms + MaxDelay_ms;
            [source1_L, source2_L] = PE_2Stims(RandStateNoiseSeq(trialnum),DUR,StimDur,DelayTimesSeq(trialnum),...
                OnsetSeq(trialnum),OffsetSeq(trialnum),whichRamp, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1),...
                LagExtSegsSeq(trialnum),DelayModsSeq(trialnum),TDT.Fs,useFile);
            
            % stim 1
            source1_R = source1_L;
		
            % remove any DCoffset
            source1_L = source1_L - mom(source1_L,1);
            source1_R = source1_R - mom(source1_R,1);
            
            %Add in the leading silent period
            source1_L =  [zeros(1,silence_len) source1_L];
            source1_R = [zeros(1,silence_len) source1_R];
            
            %Add in the trailing silent period
            source1_L =  [source1_L zeros(1,silence_len2)];
            source1_R = [source1_R zeros(1,silence_len2)];
            
            %Apply HRTF filtering
            if(XStimParams.space_flag == 1)
                if FN.HRTFfiletype(1) == 1
                    hrtf_left = mtlrch([FN.space_path FN.space],(2*LocIndSeq(1,randseq(trialnum)))-1); % !! locind -> LocIndSeq
                    hrtf_right = mtlrch([FN.space_path FN.space],2*LocIndSeq(1,randseq(trialnum)));    % !! locind -> LocIndSeq
                else
                    hrtf_left = TF1_space(LocIndSeq(1,randseq(trialnum)),:);                           % !! locind -> LocIndSeq
                    hrtf_right = TF2_space(LocIndSeq(1,randseq(trialnum)),:);                          % !! locind -> LocIndSeq
                end
                source1_L = conv(source1_L,hrtf_left);
                source1_R = conv(source1_R,hrtf_right);
            end
		
            % scale
            source1_L = source1_L - round(mean(source1_L));
            source1_R = source1_R - round(mean(source1_R));
		    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % stim 2	
            source2_R = source2_L;
         
            % remove any DCoffset
            source2_L = source2_L - mom(source2_L,1);
            source2_R = source2_R - mom(source2_R,1);
                
            %Add in the leading silent period
            source2_L =  [zeros(1,silence_len) source2_L];
            source2_R = [zeros(1,silence_len) source2_R];
            
            %Add in the trailing silent period
            source2_L =  [source2_L zeros(1,silence_len2)];
            source2_R = [source2_R zeros(1,silence_len2)];
            
            %Apply HRTF filtering
            if(XStimParams.space_flag == 1)
                if FN.HRTFfiletype(1) == 1
                    hrtf_left = mtlrch([FN.space_path FN.space],(2*LocIndSeq(2,randseq(trialnum)))-1); % !! locind -> LocIndSeq
                    hrtf_right = mtlrch([FN.space_path FN.space],2*LocIndSeq(2,randseq(trialnum)));    % !! locind -> LocIndSeq
                else
                    hrtf_left = TF1_space(LocIndSeq(2,randseq(trialnum)),:);                           % !! locind -> LocIndSeq
                    hrtf_right = TF2_space(LocIndSeq(2,randseq(trialnum)),:);                          % !! locind -> LocIndSeq
                end
                source2_L = conv(source2_L,hrtf_left);
                source2_R = conv(source2_R,hrtf_right);
            end
            
            % scale
            source2_L = source2_L - round(mean(source2_L));
            source2_R = source2_R - round(mean(source2_R));
               
            % combine stimuli
            trial_left = (source1_L + source2_L)/2;
            trial_right = (source1_R + source2_R)/2;
             
            % re-scale
            ABAval = 0.5*(mom(trial_left,2) + mom(trial_right,2));
            scalefact = TDT.scalevalue/ABAval;
            trial_left = round(scalefact*trial_left);
            trial_right = round(scalefact*trial_right);
            
            %pad with zeros
            filttrial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
            filttrial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
           
            %assignin('base','source1_R',source1_R);
            %assignin('base','source2_R',source2_R);
            if(plotAllStim | trialnum==plotStimTrial)
                PlotPEstims(source1_L, source1_R, source2_L, source2_R, TDT.Fs, ...
                    RandStateNoiseSeq(randseq(trialnum)),DelayTimesSeq(randseq(trialnum)),...
                    OnsetSeq(randseq(trialnum)),OffsetSeq(randseq(trialnum)),...
                    LagExtSegsSeq(randseq(trialnum)),DelayModsSeq(randseq(trialnum)) );
            end
		
            % save stims to disk with name of loc1
            % LocationsSeq, LocIndSeq, DelayTimesSeq, OnsetSeq, OffsetSeq, LagExtSegsSeq, DelayModsSeq
            if(exist1('H.Delayfig'));
                S232('push16',filttrial_left,length(filttrial_left));
                S232('qpop16',BUF.L1);
                %(2*LocIndSeq(2,randseq(trialnum)))-1
                fname = ['stimbuf_left_' num2str(hrtfdirmat(1,LocIndSeq(1,randseq(trialnum)))) ...
                        '_' num2str(hrtfdirmat(2,LocIndSeq(1,randseq(trialnum)))) ... 
                        '_' num2str(DelayTimesSeq(randseq(trialnum))) ...
                        '_' num2str(OnsetSeq(randseq(trialnum))) ...
                        '_' num2str(OffsetSeq(randseq(trialnum))) ...
                        '_' num2str(LagExtSegsSeq(randseq(trialnum))) ...
                        '_' num2str(DelayModsSeq(randseq(trialnum)))]; %...
                        %'_' num2str(repnum)];        % !! 2X locind -> LocIndSeq
                %disp(fname);
                evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
                        [' ''' stimuli_dir fname ''' ']   ',0);'];
                eval(evalstr);
                temp_left = dama2pc(BUF.L1);
                S232('push16',filttrial_right,length(filttrial_right));
                S232('qpop16',BUF.R1);
                fname = ['stimbuf_right_' num2str(hrtfdirmat(1,LocIndSeq(1,randseq(trialnum)))) ...
                        '_' num2str(hrtfdirmat(2,LocIndSeq(1,randseq(trialnum)))) ...
                        '_' num2str(DelayTimesSeq(randseq(trialnum))) ...
                        '_' num2str(OnsetSeq(randseq(trialnum))) ...
                        '_' num2str(OffsetSeq(randseq(trialnum))) ...
                        '_' num2str(LagExtSegsSeq(randseq(trialnum))) ...
                        '_' num2str(DelayModsSeq(randseq(trialnum)))]; %...
                        %'_' num2str(repnum)];        % !! 2X locind -> LocIndSeq
                %disp(fname);
                evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
                        [' ''' stimuli_dir fname ''' ']   ',0);'];
                eval(evalstr);
                temp_right = dama2pc(BUF.R1);
            end
            
            remtrials = numtrials - trialnum;
            set(H.remtrials,'String',num2str(remtrials));
            trialnum = trialnum + 1;
            repnum = repnum + 1;
            set(H.status,'BackgroundColor','blue');
            pause(0);
		end 	% end loop over trials
	%end         % end loop over reps
end

% if(MakePlayOption == 1)
%     return;    
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%% finished making sounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin playing sounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
set(H.status,'String','Status: Playing Stimuli');
set(H.status,'BackgroundColor','green');
set(H.remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.Delayfig') & get(H.recorddata,'Value'))
    update_dataFN;
end

%loop for reps
while (exist1('H.Delayfig') & (repnum <= XStimParams.numreps))
    % New state for each rep
    SeqRandState=sum(100*clock);
    rand('state',SeqRandState);
    % RandStateNoiseSeq= zeros(numtrials,1);
    tr=1;
    StateLow = 1; StateHi = 50; % choose state between 1 and 50
    while tr <= size(RandStateNoiseSeq,1);
        RandStateNoiseSeq(tr)=RandStates(round(StateLow + (StateHi-StateLow) * rand(1)));
        tr=tr+1;
    end
    
    %Randomize the stimuli
    if(UseRandPlaySeq)
        randseq = randperm(numtrials);  % randomized when playing
    else
        randseq = [1:numtrials];
    end
    
    trialnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    Delay_trial = [];
    OnsetMod_trial = [];
    Onset_ms_trial = [];
    LagExtSegs_ms_trial = [];
    Offset_ms_trial = [];
    RandState_trial = [];
    
    % loop for trials
    tic
    while (exist1('H.Delayfig') & (trialnum <= numtrials+1))
        
        %Check for pause by user
        if pause_check  return; end
        
        %Wait till PD1 is finished
        while S232('PD1status',1) usec_delay(1000); end
        
        if(SaveStimToDisk) % Save stimuli to disk?
            % LocationsSeq, LocIndSeq, DelayTimesSeq, OnsetSeq, OffsetSeq, LagExtSegsSeq, DelayModsSeq
            if(trialnum <= numtrials)
                fname = ['stimbuf_left_' num2str(hrtfdirmat(1,LocIndSeq(1,randseq(trialnum)))) ...
                        '_' num2str(hrtfdirmat(2,LocIndSeq(1,randseq(trialnum)))) ...
                        '_' num2str(DelayTimesSeq(randseq(trialnum))) ...
                        '_' num2str(OnsetSeq(randseq(trialnum))) ...
                        '_' num2str(OffsetSeq(randseq(trialnum))) ...
                        '_' num2str(LagExtSegsSeq(randseq(trialnum))) ...
                        '_' num2str(DelayModsSeq(randseq(trialnum)))]; %,...
                        %'_' num2str(RandStateNoiseSeq(randseq(trialnum)))   ];                                   % !! 2X locind -> LocIndSeq%
                evalstr = ['S232(''disk2dama16'',BUF.L1,'  [' ''' stimuli_dir fname ''' '] ',0);'];
                eval(evalstr);
                fname = ['stimbuf_right_' num2str(hrtfdirmat(1,LocIndSeq(1,randseq(trialnum)))) ...
                        '_' num2str(hrtfdirmat(2,LocIndSeq(1,randseq(trialnum)))) ...
                        '_' num2str(DelayTimesSeq(randseq(trialnum))) ...
                        '_' num2str(OnsetSeq(randseq(trialnum))) ...
                        '_' num2str(OffsetSeq(randseq(trialnum))) ...
                        '_' num2str(LagExtSegsSeq(randseq(trialnum))) ...
                        '_' num2str(DelayModsSeq(randseq(trialnum)))]; %,...
                        %'_' num2str(RandStateNoiseSeq(randseq(trialnum)))];                                   % !! 2X locind -> LocIndSeq%
                evalstr = ['S232(''disk2dama16'',BUF.R1,' [' ''' stimuli_dir fname ''' '] ',0);'];
                eval(evalstr);
            end
        else
            % ////////////////////////////////////////////
             %Make the stimuli on the fly
            if(trialnum <= numtrials)
                %Make the stimuli
                useFile=0;
%                 if get(H.Delay_stim_type,'Value') == 2  % Generate
%                     useFile=1;
%                 end
                whichRamp=get(H.DelayOnOff_mode,'Value'); % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
                
                %StimDur = DUR + MaxLag_ms + MaxDelay_ms;
                      
                [source1_L, source2_L] = PE_2Stims(RandStateNoiseSeq(randseq(trialnum)),DUR,StimDur,DelayTimesSeq(randseq(trialnum)),...
                    OnsetSeq(randseq(trialnum)),OffsetSeq(randseq(trialnum)),whichRamp, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1),...
                    LagExtSegsSeq(randseq(trialnum)),DelayModsSeq(randseq(trialnum)),TDT.Fs,useFile);
                
                % stim 1
                source1_R = source1_L;
			
                % remove any DCoffset
                source1_L = source1_L - mom(source1_L,1);
                source1_R = source1_R - mom(source1_R,1);
                
                %Add in the leading silent period
                source1_L =  [zeros(1,silence_len) source1_L];
                source1_R = [zeros(1,silence_len) source1_R];
                
                %Add in the trailing silent period
                source1_L =  [source1_L zeros(1,silence_len2)];
                source1_R = [source1_R zeros(1,silence_len2)];
                
                %Apply HRTF filtering
                if(XStimParams.space_flag == 1)
                    if FN.HRTFfiletype(1) == 1
                        hrtf_left = mtlrch([FN.space_path FN.space],(2*LocIndSeq(1,randseq(trialnum)))-1); % !! locind -> LocIndSeq
                        hrtf_right = mtlrch([FN.space_path FN.space],2*LocIndSeq(1,randseq(trialnum)));    % !! locind -> LocIndSeq
                    else
                        hrtf_left = TF1_space(LocIndSeq(1,randseq(trialnum)),:);                           % !! locind -> LocIndSeq
                        hrtf_right = TF2_space(LocIndSeq(1,randseq(trialnum)),:);                          % !! locind -> LocIndSeq
                    end
                    source1_L = conv(source1_L,hrtf_left);
                    source1_R = conv(source1_R,hrtf_right);
                end
			
                % scale
                source1_L = source1_L - round(mean(source1_L));
                source1_R = source1_R - round(mean(source1_R));
			    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % stim 2	
                source2_R = source2_L;
             
                % remove any DCoffset
                source2_L = source2_L - mom(source2_L,1);
                source2_R = source2_R - mom(source2_R,1);
                    
                %Add in the leading silent period
                source2_L =  [zeros(1,silence_len) source2_L];
                source2_R = [zeros(1,silence_len) source2_R];
                
                %Add in the trailing silent period
                source2_L =  [source2_L zeros(1,silence_len2)];
                source2_R = [source2_R zeros(1,silence_len2)];
                
                %Apply HRTF filtering
                if(XStimParams.space_flag == 1)
                    if FN.HRTFfiletype(1) == 1
                        hrtf_left = mtlrch([FN.space_path FN.space],(2*LocIndSeq(2,randseq(trialnum)))-1); % !! locind -> LocIndSeq
                        hrtf_right = mtlrch([FN.space_path FN.space],2*LocIndSeq(2,randseq(trialnum)));    % !! locind -> LocIndSeq
                    else
                        hrtf_left = TF1_space(LocIndSeq(2,randseq(trialnum)),:);                           % !! locind -> LocIndSeq
                        hrtf_right = TF2_space(LocIndSeq(2,randseq(trialnum)),:);                          % !! locind -> LocIndSeq
                    end
                    source2_L = conv(source2_L,hrtf_left);
                    source2_R = conv(source2_R,hrtf_right);
                end
                
                % scale
                source2_L = source2_L - round(mean(source2_L));
                source2_R = source2_R - round(mean(source2_R));
                   
                % combine stimuli
                trial_left = (source1_L + source2_L)/2;
                trial_right = (source1_R + source2_R)/2;
                
                %///////////////////////////////////////
                % Source 1 alone for testing !!!!!!!!!!!!
%                 if(0)
%                     trial_left = source1_L;
%                     trial_right = source1_R;
%                 end
                 
                % re-scale
                ABAval = 0.5*(mom(trial_left,2) + mom(trial_right,2));
                scalefact = TDT.scalevalue/ABAval;
                trial_left = round(scalefact*trial_left);
                trial_right = round(scalefact*trial_right);
                
%                 figure;
%                 plot(trial_left);
%                 return;
                
                %pad with zeros
                filttrial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
                filttrial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
                
                S232('push16',filttrial_left,length(filttrial_left));
                S232('qpop16',BUF.L1);
                S232('push16',filttrial_right,length(filttrial_right));
                S232('qpop16',BUF.R1);
                
                %evalstr = ['S232(''disk2dama16'',BUF.R1,' [' ''' stimuli_dir fname ''' '] ',0);'];
                %eval(evalstr);
                
                %assignin('base','source1_R',source1_R);
                %assignin('base','source2_R',source2_R);
                if(plotAllStim | trialnum==plotStimTrial)
                    PlotPEstims(source1_L, source1_R, source2_L, source2_R, TDT.Fs, ...
                        RandStateNoiseSeq(randseq(trialnum)),DelayTimesSeq(randseq(trialnum)),...
                        OnsetSeq(randseq(trialnum)),OffsetSeq(randseq(trialnum)),...
                        LagExtSegsSeq(randseq(trialnum)),DelayModsSeq(randseq(trialnum)) );
                end
            end
        end
        
        %Wait till PD1 is finished
        while S232('PD1status',1) usec_delay(1000); end
        S232('PD1stop',1);
        
        %Stop the m110 and get spikes
        if(trialnum > 1)							% first trial just for loading
            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, 3*XStimParams.curr_stimdur); 			% Take XStimParams.curr_stimdur spikes max +50 ms!!!!
            ind = find(spikes ~= 0); 						% Get clock events that are spikes
            spikes = spikes(ind);
            ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
            if(exist1('H.Delayfig') & ~isempty(spikes)) 
                spikes = spikes(ind);
                spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
                EL_trial = [EL_trial;LocationsSeq(1,randseq(trialnum-1))* ones(size(spikes))];
                AZ_trial = [AZ_trial;LocationsSeq(2,randseq(trialnum-1))* ones(size(spikes))];
                repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
                Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
                
                Delay_trial = [Delay_trial;DelayTimesSeq(randseq(trialnum-1),1)* ones(size(spikes))];
                Onset_ms_trial = [Onset_ms_trial;OnsetSeq(randseq(trialnum-1),1)* ones(size(spikes))];
                Offset_ms_trial = [Offset_ms_trial;OffsetSeq(randseq(trialnum-1),1)* ones(size(spikes))];
                LagExtSegs_ms_trial = [LagExtSegs_ms_trial;LagExtSegsSeq(randseq(trialnum-1),1)* ones(size(spikes))];
                OnsetMod_trial = [OnsetMod_trial;DelayModsSeq(randseq(trialnum-1),1)* ones(size(spikes))];
                RandState_trial = [RandState_trial;RandStateNoiseSeq(randseq(trialnum-1),1)* ones(size(spikes))];
            end
        end
        
        %Check for pause by user
        if pause_check  return; end

        if(exist1('H.Delayfig') & (trialnum <= numtrials))
            
            S232('seqplay',BUF.playspec1);
            S232('PD1arm',1);
            
            %Send trigger
            %Set up MII
            m100x( C_.INIT );
            m110dx( C_.INIT );
            m110dx( C_.CLOCK, mii_us_per_sample);
            m110dx( C_.MODE, M_.PST );
            
            if (trialnum <= numtrials)
                while toc < ISI/1000     end
                %Start clock
                m110dx( C_.START);
                %Send pulse: PD1 GO!
                m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
                tic
            end
            
        end
        
        if(trialnum > 1)
            finalspikematrix(randseq(trialnum-1)) = ...
                finalspikematrix(randseq(trialnum-1)) + ...
                length(spikes);
        end
        
        remtrials = numtrials - trialnum +1;
        set(H.remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end %end loop over trials
    
    if(1)
        if(repnum<=20)
            %Plot Spike Rate Data
            if(findobj('Tag','Raster'))
                %close('Raster');
                figure(H.delay_raster);
                title(['Rep# ' num2str(repnum)]);
                %ylabel('trial'); xlabel('time (ms)'); title(['Rep# ' num2str(repnum)]);
                %ylim([0 numtrials+1]);
                %set(H.delay_raster_axes, 'YLim', [0 numtrials]);
            else
				H.delay_raster = figure('Units','pixels',...
                  'Position',[870 0 400 920],...
                  'Tag', 'Raster',...
                  'Name','Raster',...
                  'NumberTitle','off',...
                  'MenuBar', 'none',...
                  'Color',[0 0 0]);
                H.delay_raster_axes = axes;
                figure(H.delay_raster);
                ylabel('trial'); xlabel('time (ms)'); title(['Rep# ' num2str(repnum)]);
                ylim([0 numtrials+1]);
                xmax=XStimParams.silence_lead +DUR +50 +MaxLag_ms; %MaxDelay_ms
                xlim([XStimParams.silence_lead xmax])
                set(H.delay_raster_axes, 'YLim', [0 numtrials+0.5]);
                set(H.delay_raster_axes, 'GridLineStyle', ':');
                set(H.delay_raster_axes, 'XGrid', 'on');
                set(H.delay_raster_axes, 'Color', [0 0 0]);
                set(H.delay_raster_axes, 'XColor', [0 0 1]);
                set(H.delay_raster_axes, 'YColor', [0 0 1]);
            end
            plotmatrix = [Nspikes spikes_trial repnum_trial EL_trial AZ_trial Delay_trial  Onset_ms_trial Offset_ms_trial LagExtSegs_ms_trial OnsetMod_trial];
            %assignin('base','plotmatrix',plotmatrix);
            Delay_raster(XStimParams, plotmatrix, randseq, XStimParams.numreps);
        end
        %Record Data
        if(exist1('H.Delayfig') & get(H.recorddata,'Value') )
            tempseq{repnum} = randseq;
            datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial Delay_trial  Onset_ms_trial Offset_ms_trial LagExtSegs_ms_trial OnsetMod_trial RandState_trial]];
            record_data3(XStimParams,datamatrix, tempseq);
            % output: el az delay_ms onset_ms offset_ms onset_mod rep num_spikes spikes[]
            % onset_mod: 1=normal, 2=uncorr, 3=corr without onset
        end
    end
    remreps = XStimParams.numreps - repnum;
    set(H.remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

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

% ////////////////////////////////////
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
function [stim] = PEstim(state,SndDur,Fs,UseFiles,flip)
    % PEstim: Creates or loads a BroadBand Noise (2-11 kHz) with
    % specified: rand state (state) and sound duration.
    % passing a state of 0 loads stim 'randTest' and verifies rand #
    % generator
    
    %filePath = 'C:\Brian\PEStim\';
    filePath = 'E:\Brian\PEStim\';
    SaveFile=0; % Set to 1 manually when making new noises
    
    if nargin < 4 UseFiles=0; end
    if nargin < 3 Fs=30000; end
    if nargin < 2 SndDur=300; end
    
    if(SndDur>300)
        warning('Noise files are only 300 ms long....');    
    end
    
    % make sure rand # generator hasn't changed, use files if it has...
    if state == 0 % option for saving RandTest
        state=999; % state for testing
        RandTest = GetNoise_BB(state,30000,1);
        %round(std(RandTest)*100000))= 41824
        state=0;
        if(SaveFile)
            save([filePath 'randGenTest'], 'RandTest');
            return;
        end
    end
    if UseFiles | state == 0  % load saved file corresponding with state number
        if state ~= 0
            fileName = ['PEStim_' num2str(state)];
            load([filePath fileName], 'stim');
            % truncate to SndDur
            pnts=(SndDur/1000)*Fs+1;
            stim_pnts=size(stim,2);
            stim(pnts:stim_pnts)=[];
            
        else % test rand generator
           fileName = 'randGenTest';
           load([filePath fileName], 'RandTest');
           %disp(round(std(RandTest)*100000));
           if(round(std(RandTest)*100000)== 41824)
               stim=1;
           else
               stim=0;
           end
           return;
        end
        
        if(size(stim,1) > 1)
            warning 'stim loaded in wrong dimension';
            stim = stim';
        end
        
    else    % make stimulus using state    
        stim = GetNoise_BB(state,Fs,300);
        if(SaveFile)
            fileName = ['PEStim_' num2str(state)];
            save([filePath fileName], 'stim');
            return;
        end
        if(flip) % uncorrelated
            stim = fliplr(stim); % transpose
        end
        % truncate to SndDur
        pnts=(SndDur/1000)*Fs+1;
        stim_pnts=size(stim,2);
        stim(pnts:stim_pnts)=[];
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


% //////////////////////////////////////////////////   
function [stim1, stim2] = PE_2Stims(state,SndDur,StimDur,delay,onset,offset,whichRamp,onDef,offDef,LagExt,Mod,Fs,useFile)
    % PE_2Stims: Loads/Generates, Delays, ramps, and buffers stim
    
    % Stim 1
    % Lead of Lag?
    if(delay<=0)   % is lead
        Delay_1=0;
        LagExt_1=0;
    else            % is lag
        Delay_1=delay;
        LagExt_1=LagExt;
    end
   
    SilentOnset=0;
    SilentOffset=0;
    % Stim duration mods, modify only the lagging sound
    if((Mod==3 | Mod==5) & delay < 0) % silence lead segment
        SilentOnset=1;
        Delay_1=delay;  % must send in delay if silencing the lead segment
    else
        SilentOnset=0;
    end
    if(Mod==4 & delay > 0) % silence lag segment
        SilentOffset=1;
        Delay_1=delay; % must send in delay if silencing its laggin seggment
        LagExt_1=0;
    else
        SilentOffset=0;
    end
    flip=0; % never flip stim1 (only stim2)
    
    % Get sound with appropriate duration
    stim1 = PEstim(state,SndDur+LagExt_1,Fs,useFile,flip);
    pnts_1=size(stim1,2);

    if SilentOnset % 1 or 0  % silence lead segment
        Dpnts=round((abs(delay)/1000)*Fs-1);
        stim1(1:Dpnts)=[];
    end

    if SilentOffset % 1 or 0  % silence lag segment
        Offpnts=round((abs(delay)/1000)*Fs-1);
        stim1(pnts_1-Offpnts:pnts_1)=[];
        %Delay_1=Delay_1*2;
    end
    pnts_1=size(stim1,2);
    
    % ABS delay (stim #2 is delayed when delay is negative)
    Dpnts=round((abs(Delay_1)/1000)*Fs-1);

    % total stimulus duration
    StimPnts=round((StimDur/1000)*Fs);
    if(StimPnts < Dpnts+pnts_1)
        StimPnts = Dpnts+pnts_1;
        warning('The sound and its delay are longer than stimulus duration');
    end

    %ramp stim   whichRamp
    switch whichRamp % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
        case 1
            stim1 = ramp_sound_onoff(stim1,Fs,onset,offset);
        case 2
            stim1 = ramp_sound_onoff(stim1,Fs,onset,offset);
        case 3
            stim1 = ramp_sound_onoff(stim1,Fs,onDef,offDef);
        otherwise
            disp('Onset Mode set to 1 - source 1 alone uses random onsets and offsets');
            stim1 = ramp_sound_onoff(stim1,Fs,onset,offset);
    end
    
    % buffer with ABS delay
    if(Dpnts>0)
        Silence=zeros(1,Dpnts+pnts_1); % stimulus duration
        Silence(size(Silence,2)-size(stim1,2)+1:size(Silence,2))=stim1;
        stim1 = Silence;
    end
    
    % buffer with ABS delay
    if(Dpnts>0)
        Silence=zeros(1,Dpnts+pnts_1); % stimulus duration
        Silence(size(Silence,2)-size(stim1,2)+1:size(Silence,2))=stim1;
        stim1 = Silence;
    end
    
    
    % Stim 2
    % Lead of Lag?
    if(delay>0)   % is lead
        Delay_2=0;
        LagExt_2=0;
    else            % is lag
        Delay_2=delay;
        LagExt_2=LagExt;
    end
    
    SilentOnset=0;
    SilentOffset=0;
    % Stim duration mods, modify only the lagging sound
    if((Mod==3 | Mod==5) & delay > 0) % silence lead segment
        SilentOnset=1;
        Delay_2=delay;  % must send in delay if silencing the lead segment
    else
        SilentOnset=0;
    end
    if(Mod==4 & delay < 0) % silence lag segment
        SilentOffset=1;
        Delay_2=delay; % must send in delay if silencing its laggin seggment
        LagExt_2=0;
    else
        SilentOffset=0;
    end
    % flip noise if uncorrelated (it a long (300ms noise).
    flip=0;
    if(Mod==2 | Mod==5)
        flip=1;
    end
    
     % Get sound with appropriate duration
    stim2 = PEstim(state,SndDur+LagExt_2,Fs,useFile,flip);
    pnts_2=size(stim2,2);

    if SilentOnset % 1 or 0  % silence lead segment
        Dpnts=round((abs(delay)/1000)*Fs-1);
        %disp([Dpnts length(stim2)]);
        stim2(1:Dpnts)=[];
    end

    if SilentOffset % 1 or 0  % silence lag segment
        Offpnts=round((abs(delay)/1000)*Fs-1);
        stim2(pnts_2-Offpnts:pnts_2)=[];
        %Delay_2=Delay_2*2;
    end
    pnts_2=size(stim2,2);
    
    % ABS delay (stim #2 is delayed when delay is negative)
    Dpnts=round((abs(Delay_2)/1000)*Fs-1);

    % total stimulus duration
    StimPnts=round((StimDur/1000)*Fs);
    if(StimPnts < Dpnts+pnts_2)
        StimPnts = Dpnts+pnts_2;
        warning('The sound and its delay are longer than stimulus duration');
    end

    %ramp stim   whichRamp
    switch whichRamp % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
        case 1
            stim2 = ramp_sound_onoff(stim2,Fs,onDef,offDef);
        case 2
            stim2 = ramp_sound_onoff(stim2,Fs,onset,offset);
        case 3
            stim2 = ramp_sound_onoff(stim2,Fs,onset,offset);
        otherwise
            disp('Onset Mode set to 1 - source 1 alone uses random onsets and offsets');
            stim2 = ramp_sound_onoff(stim2,Fs,onDef,offDef);
    end
    
    % buffer with ABS delay
    if(Dpnts>0)
        Silence=zeros(1,Dpnts+pnts_2); % stimulus duration
        Silence(size(Silence,2)-size(stim2,2)+1:size(Silence,2))=stim2;
        stim2 = Silence;
    end
    
   
    % buffer to stimulus' duration
    Silence=zeros(1,StimPnts);
    Silence(1:size(stim1,2))=stim1;
    stim1 = Silence;
    
    Silence=zeros(1,StimPnts);
    Silence(1:size(stim2,2))=stim2;
    stim2 = Silence;
    
    if(0)
        timeaxis = 0 : 1/Fs : (size(stim1,2)*(1/Fs))-(1/Fs);
        if(findobj('Tag','Stimuli')) close('Stimuli');  end
        Stim.fig = figure('Units','pixels',...
          'Position',[50 500 500 500],...
          'Tag', 'Stimuli',...
          'Name','Stimuli',...
          'NumberTitle','off',...
          'Color',[1 1 1]);
        hold on;
        Stim.L1 = plot(timeaxis, stim1 + 0.9,...
            'Tag','Env_ellipses',...
            'Marker', 'none',...
            'LineStyle','-',...
            'Color', [1 0 0],...
            'LineWidth', 0.5);
        Stim.L2 = plot(timeaxis, stim2 - 0.9,...
            'Tag','Env_ellipses',...
            'Marker', 'none',...
            'LineStyle','-',...
            'Color', [0 0 1],...
            'LineWidth', 0.5);
        hold off;
        title('red = 1 (target), blue = 2 (masker)');
    end
    
   
% //////////////////////////////////////////
function PlotPEstims(source1_L, source1_R, source2_L, source2_R, Fs ,p1, p2, p3, p4, p5, p6)
        timeaxis = 0 : 1/Fs : (size(source1_L,2)*(1/Fs))-(1/Fs);
        %disp(size(timeaxis,2));
        %disp(size(source1_L,2));
        %return;
        if(findobj('Tag','Stimuli'))
            close('Stimuli');
        end
        Stim.fig = figure('Units','pixels',...
          'Position',[50 500 500 500],...
          'Tag', 'Stimuli',...
          'Name','Stimuli',...
          'NumberTitle','off',...
          'Color',[1 1 1]);
        hold on;
        %subplot(2,1,1);
        Stim.R = plot(timeaxis, source1_R + 0.9,...
            'Tag','Env_ellipses',...
            'Marker', 'none',...
            'LineStyle','-',...
            'Color', [1 0.5 0],...
            'LineWidth', 0.5);
        Stim.L = plot(timeaxis, source1_L + 0.9,...
            'Tag','Env_ellipses',...
            'Marker', 'none',...
            'LineStyle','-',...
            'Color', [1 0 0],...
            'LineWidth', 0.5);
        %subplot(2,1,2);
        Stim.R = plot(timeaxis, source2_R - 0.9,...
            'Tag','Env_ellipses',...
            'Marker', 'none',...
            'LineStyle','-',...
            'Color', [0 0.5 1],...
            'LineWidth', 0.5);
        Stim.L = plot(timeaxis, source2_L - 0.9,...
            'Tag','Env_ellipses',...
            'Marker', 'none',...
            'LineStyle','-',...
            'Color', [0 0 1],...
            'LineWidth', 0.5);
        hold off;
        %set(gca, 'Title', 'red=src1, blue=src2');
        title('Source: red=1 (target), blue=2 (masker)');
        %text(0.3,0.3,'[one two]');

        str = strvcat(num2str(p1), num2str(p2), num2str(p3), num2str(p4), num2str(p5), num2str(p6));
        text(size(source1_L,2)*(1/Fs)+0.01, 1, str);
