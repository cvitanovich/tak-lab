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
% XStimParams added to XStimParams_func  Run to add!!!

% to do:
% add LAG silence back in???

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
    %disp(DelayTimes);
    MaxDelay_pnts=round(TDT.Fs*(max1(DelayTimes)/1000));
    
    % Get rid of NaN onset ramps
    ramp_timeS_on=zeros(1,1);
    mp=1;
    for mt = 1:size(XStimParams.ramp_timeS_on,2)
        if(~isnan(XStimParams.ramp_timeS_on(mt)))
            ramp_timeS_on(mp)=XStimParams.ramp_timeS_on(mt);
            mp=mp+1;
        end
    end
    MaxOnset=max1(ramp_timeS_on);
    %disp(ramp_timeS_on);
    
    % Get rid of NaN offset ramps
    ramp_timeS_off=zeros(1,1);
    mp=1;
    for mt = 1:size(XStimParams.ramp_timeS_off,2)
        if(~isnan(XStimParams.ramp_timeS_off(mt)))
            ramp_timeS_off(mp)=XStimParams.ramp_timeS_off(mt);
            mp=mp+1;
        end
    end
    MaxOffset=max1(ramp_timeS_off); % ramp_timeS_off ???????
    %disp(ramp_timeS_off);
    
    lag_seg_ext=zeros(1,1);
    mp=1;
    for mt = 1:size(XStimParams.lag_seg_ext,2)
        if(~isnan(XStimParams.lag_seg_ext(mt)))
            lag_seg_ext(mp)=XStimParams.lag_seg_ext(mt);
            mp=mp+1;
        end
    end
    MaxDUR=max1(lag_seg_ext);  
    %disp(lag_seg_ext);
  
else
    SilentDelay=XStimParams.SilentDelay;
    if(SilentDelay & XStimParams.DelayTimes(1) >= XStimParams.curr_stimdur)
        XStimParams.DelayTimes(1)=XStimParams.curr_stimdur-10;
        set(H.delaytime,'String', num2str(XStimParams.DelayTimes(1)));
        warndlg('delay set to stim. dur. minus 10 ms', 'warning!');
        waitforbuttonpress;
    end
    Delay_pnts=round(TDT.Fs*(XStimParams.DelayTimes(1)/1000));
    uncorrel=XStimParams.uncorrel;
    SilentLag=XStimParams.SilentLag;
    MaxDelay_pnts=Delay_pnts;
    
    ramp_time_on = XStimParams.ramp_timeS_on(1); %ms
    ramp_time_off = XStimParams.ramp_timeS_off(1); %ms
    MaxDUR = XStimParams.lag_seg_ext(1);
end

if(findobj('Tag','Raster'))
    close('Raster');
end

rand('state',sum(100*clock));  % Fresh stimulus
 
stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);

%check if filt files assigned
if exist([FN.ephone_path FN.ephone2]) ~= 2
    ephonefilediagbox;
end
if XStimParams.space_flag |  XStimParams.ABLalone_flag
    while exist1([FN.space_path FN.space]) ~=2
        [FN.space,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select Fully-cued HRTF File');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space]);
        end
        FN.HRTFfiletype(1) = testHRTFfiletype(FN.space_path, FN.space);
    end
    if XStimParams.space_flag
        disp('This is a FULLY-CUED SPACE test')
    else
        disp('This is an ABL-alone SPACE test')
    end
elseif XStimParams.ildalone_flag
    while exist1([FN.ILA_path FN.ildalone]) ~=2
        [FN.ildalone,FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.ildalonefile,'String',[FN.ILA_path FN.ildalone]);
        end      
        FN.HRTFfiletype(2) = testHRTFfiletype(FN.ILA_path, FN.ildalone);
    end
    disp('This is an ILD-ALONE test')
elseif XStimParams.itdalone_flag
    while exist1([FN.ILA_path FN.ildalone]) ~=2
        [FN.ildalone,FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.ildalonefile,'String',[FN.ILA_path FN.ildalone]);
        end
        FN.HRTFfiletype(2) = testHRTFfiletype(FN.ILA_path, FN.ildalone);
    end
    while exist1([FN.ITA_path FN.itdalone]) ~=2
        [FN.itdalone,FN.ITA_path] = uigetfile([FN.ITA_path '*.*'],'Select ITD-alone HRTF File');
        if(FN.ITA_path ~= 0)
            set(H.itdalonefile,'String',[FN.ITA_path FN.itdalone]);
        end
        FN.HRTFfiletype(3) = testHRTFfiletype(FN.ITA_path, FN.itdalone);
    end
    disp('This is an ITD-ALONE test')
end

%Put parameters into XStimParams
XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.reset_flag = 0;

if XStimParams.ABLalone_flag
    % load filter coeffs & other params for ABLalone test
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

%Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
silence_len = (XStimParams.silence_lead * round(TDT.Fs/1000));
%Add a piece of silence after stimulus 
silence_len2 = (XStimParams.silence_trail * round(TDT.Fs/1000));

%Make Stimulus buffers
DUR = XStimParams.curr_stimdur + MaxDUR; % add MaxDUR to account for longer lag
if XStimParams.ildalone_flag 					%ILDAlone Test
    S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen + MaxDelay_pnts);
    S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen + MaxDelay_pnts);
    S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen + MaxDelay_pnts);
    S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen + MaxDelay_pnts);
elseif XStimParams.itdalone_flag 				%ITDAlone Test
    S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2 + MaxDelay_pnts);
    S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2 + MaxDelay_pnts);
    S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2 + MaxDelay_pnts);
    S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2 + MaxDelay_pnts);
elseif XStimParams.space_flag | XStimParams.ABLalone_flag	    %fully-cued Test
    S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen + MaxDelay_pnts);
    S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen + MaxDelay_pnts);
    S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen + MaxDelay_pnts);
    S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen + MaxDelay_pnts);
end   

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
if XStimParams.ildalone_flag 					%ILDAlone Test
    S232('PD1npts',1,(silence_len2 + silence_len + DUR*(round(TDT.Fs/1000))) + TDT.itdfiltlen  + TDT.ephonefiltlen + TDT.hrtffiltlen + MaxDelay_pnts);
elseif XStimParams.itdalone_flag 				%ITDAlone Test
    S232('PD1npts',1,(silence_len2 + silence_len + DUR*(round(TDT.Fs/1000)))  + TDT.ephonefiltlen + (2 * TDT.hrtffiltlen) + MaxDelay_pnts);
elseif XStimParams.space_flag |  XStimParams.ABLalone_flag			% fully-cued or ABL-alone Test
    S232('PD1npts',1,(silence_len2 + silence_len + DUR*(round(TDT.Fs/1000))) + TDT.ephonefiltlen + TDT.hrtffiltlen + MaxDelay_pnts);
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

if XStimParams.space_flag | XStimParams.ABLalone_flag
    if FN.HRTFfiletype(1) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space]));
    elseif FN.HRTFfiletype(1) == 2
        hrtfdirmat = dir_space;
    end
elseif XStimParams.ildalone_flag
    if FN.HRTFfiletype(2) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone]));
    elseif FN.HRTFfiletype(2) == 2
        hrtfdirmat = dir_ILA;
    end
elseif XStimParams.itdalone_flag
    if FN.HRTFfiletype(3) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.ITA_path FN.itdalone]));
    elseif FN.HRTFfiletype(3) == 2
        hrtfdirmat = dir_ITA;
    end
end
%%%%%%%%%%%

% find locations for 2 or 3 sources
clear locind locations
numlocs = 0;

offset_el = XStimParams.offset_el(2);%/2;
offset_az = XStimParams.offset_az(2);%/2;
% offset_el3 = XStimParams.offset_el(3)/2;
% offset_az3 = XStimParams.offset_az(3)/2;
% if ~XStimParams.ThreeStims
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
% else
%     for locnum = 1:size(XStimParams.locations,2)
%         loc1 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,locnum) + offset_el) &...
%             hrtfdirmat(2,:) == (XStimParams.locations(2,locnum)) + offset_az));
%         loc2 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,locnum) - offset_el) &...
%             hrtfdirmat(2,:) == (XStimParams.locations(2,locnum)) - offset_az));
%         loc3 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,locnum) + offset_el3) &...
%             hrtfdirmat(2,:) == (XStimParams.locations(2,locnum)) + offset_az3));
%         if ~isempty(loc1) & ~isempty(loc2) & ~isempty(loc3)
%             numlocs = numlocs+1;
%             locind(:,numlocs) = [loc1; loc2; loc3];
%             locations(:,numlocs) = XStimParams.locations(:,locnum);
%         end
%     end
%     
%end

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
    %  = [1 2 3]; % 1=normal, 2=uncorr, 3=corr without onset
    for tr = 1:numtrials
        if(uncorrel & ~SilentDelay & ~SilentLag)
            DelayModsSeq(tr)=2;
        else
            if(~uncorrel & SilentDelay & ~SilentLag)
                DelayModsSeq(tr)=3;
            else
                if(~uncorrel & ~SilentDelay & SilentLag)
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
    
else % random parameters
    % getting these here because they don't seem to be updated properly...
    XStimParams.uncorrel = get(H.uncorrel,'Value');
    uncorrel=XStimParams.uncorrel;
    XStimParams.SilentDelay = get(H.silentdelay,'Value');
    SilentDelay=XStimParams.SilentDelay;
    XStimParams.SilentLag = get(H.silentlag,'Value');
    SilentLag=XStimParams.SilentLag;
    
    if(uncorrel & SilentDelay & SilentLag)
        XStimParams.DelayMods = [1 2 3 4]; 
        % 1=normal, 2=uncorr, 3=corr no lead segment, 4=corr no lag segment
    else
        if(~uncorrel & SilentDelay & ~SilentLag)
            XStimParams.DelayMods = [1 3]; % 1=normal, 3=corr without onset
        else
            if(uncorrel & ~SilentDelay & ~SilentLag)
                XStimParams.DelayMods = [1 2]; % 1=normal
            else
                if(~uncorrel & ~SilentDelay & SilentLag)
                    XStimParams.DelayMods = [1 4]; % 1=normal
                else
                   if(~uncorrel & ~SilentDelay & ~SilentLag)
                        XStimParams.DelayMods = [1]; % 1=normal
                   else
                       warning('no valid stimulus combination selected');
                       XStimParams.DelayMods = [1]; % 1=normal
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
    
    % LocationsSeq, LocIndSeq, DelayTimesSeq, OnsetSeq, OffsetSeq, LagExtSegsSeq, DelayModsSeq
     
    % Eliminate trials where delay exceeds stim. duration - on/offset but is silenced
    OrigNumtrials=numtrials;
    numtrials=1;
    tr=1;
    while tr <= size(DelayModsSeq,1);
        if(DelayModsSeq(tr)==3 & DelayTimesSeq(tr) > XStimParams.curr_stimdur-(MaxOnset+MaxOffset) )
            LocationsSeq(:,tr)=[,];
            LocIndSeq(:,tr)=[,];
            DelayTimesSeq(tr)=[];
            OnsetSeq(tr)=[];
            OffsetSeq(tr)=[];
            LagExtSegsSeq(tr)=[];
            DelayModsSeq(tr)=[];
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
    
    % Eliminate trials where DelayModsSeq == 4 (ext. lag) but the lag is silenced
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
        else
            tr=tr+1;
            numtrials=numtrials+1;
        end
    end
    numtrials=numtrials-1;

    disp(['trial ',...
    'elevation ',...
    'azimuth ',...
    'delay ',...
    'onset ',...
    'offset ',...
    'LagSegExt ',...
    'mod (1=async 2=uncorr 3=no lead seg 4=no lag seg)'...
    ]);

    SeqOrder = zeros(7,numtrials);
    for d = numtrials:-1:1
        SeqOrder(1,d) = d;
        SeqOrder(2,d) = LocationsSeq(1,d);
        SeqOrder(3,d) = LocationsSeq(1,d);
        SeqOrder(4,d) = DelayTimesSeq(d);
        SeqOrder(5,d) = OnsetSeq(d);
        SeqOrder(6,d) = OffsetSeq(d);
        SeqOrder(7,d) = LagExtSegsSeq(d);
        SeqOrder(8,d) = DelayModsSeq(d);
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
            ]);
    end
end

finalspikematrix = zeros(1,numtrials);

%Randomize the stimuli
randseq = randperm(numtrials);
trialnum = 1;
while (exist1('H.Delayfig') & (trialnum <= numtrials))
    set(H.status,'BackgroundColor','red');
    %Check for pause by user
    if pause_check  return; end
    
   %Make the first stimulus
    if get(H.stim_type,'Value') ~= 9
        source1_L = get_stim(XStimParams.freq(1)); % adding 50 ms in get_stim !!!
    else
        source1_L = get_stim(FN.stim_path,FN.stim);
    end
        
    source1_R = source1_L;
    
    % copy fist source in case not generated
    if(~uncorrel)
        source2_L = source1_L;
        source2_R = source1_R;
    end
    
    origPnts=size(source1_L,2);
    
    % delete 50-(50-LagExtSegsSeq(trialnum)) of stim since I added 50 ms in get_stim
    if(LagExtSegsSeq(trialnum)~=0)
        ExtDel_pnts=round(TDT.Fs*((50-LagExtSegsSeq(trialnum)))/1000); % added 50 ms in get_stim
    else
        ExtDel_pnts=round(TDT.Fs*(50)/1000); % added 50 ms in get_stim
    end
    source1_L(size(source1_L,2)-ExtDel_pnts : size(source1_L,2))=[];
    source1_R(size(source1_R,2)-ExtDel_pnts : size(source1_R,2))=[];
    ExtDel_pnts=ExtDel_pnts+1;
    
    % Silence lag if DelayModsSeq = 4
    SegDel_pnts=0;
    if(DelayModsSeq(trialnum)==4)
        SegDel_pnts=round(TDT.Fs*(DelayTimesSeq(trialnum))/1000); % added 50 ms in get_stim
        source1_L(size(source1_L,2)-SegDel_pnts : size(source1_L,2))=[];
        source1_R(size(source1_R,2)-SegDel_pnts : size(source1_R,2))=[];
        SegDel_pnts=SegDel_pnts+1;
        %SegDel_pnts=abs(SegDel_pnts);
    end
    
    % modulate stim1
    % DUR IS PROBABLY NO LONGER VALID HERE!!!  
    if ~strcmp(XStimParams.mod_type,'None')
        Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
        source1_R = source1_R .* Envelope(:)';
        source1_L = source1_L .* Envelope(:)';
    end
    
     % set onset type here so that it doesn't have to be random
    if(XStimParams.randOnsetPerms)
        if(DelayModsSeq(trialnum)==2)
            uncorrel=1;
        else
            uncorrel=0;
        end
    end

    % get delay time, get correl/uncorrel, get onset type
    % will buffer to MaxDelay_Pnts
    if(XStimParams.randOnsetPerms)
        Delay_pnts=round(TDT.Fs*(DelayTimesSeq(trialnum))/1000);
        if(Delay_pnts > 0)
            MaxDelayBufferPnts = MaxDelay_pnts - Delay_pnts;
        else
            MaxDelayBufferPnts = MaxDelay_pnts + Delay_pnts;
        end
    end
    
    % set onset type here so that it doesn't have to be random
    if(XStimParams.randOnsetPerms)
        if(DelayModsSeq(trialnum)==3)
            SilentDelay=1;
        else
            SilentDelay=0;
        end
    end
    
    if(Delay_pnts<0 & SilentDelay) % will zero time of delay in first stim. SEE BELOW
        origPnts=size(source1_L,2);
        
        source1_L_temp = source1_L;
        source1_R_temp = source1_R;
        
        % delete delay
        %disp(size(source1_L_temp,2));
        source1_L = source1_L_temp(-Delay_pnts:size(source1_L_temp,2));
        source1_R = source1_R_temp(-Delay_pnts:size(source1_R_temp,2));
        %disp(size(source1_L,2));
        
        %if get(H.stim_type,'Value') ~= 9      %Ramp the stimuli
            %ramp_time_on = 10; %ms
            %ramp_time_off = 2.5; %ms
            switch get(H.DelayOnOff_mode,'Value') % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
                case 1
                    [source1_L] = ramp_sound_onoff(source1_L,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                    [source1_R] = ramp_sound_onoff(source1_R,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                case 2
                    [source1_L] = ramp_sound_onoff(source1_L,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                    [source1_R] = ramp_sound_onoff(source1_R,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                case 3
                    [source1_L] = ramp_sound_onoff(source1_L,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
                    [source1_R] = ramp_sound_onoff(source1_R,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
                otherwise
                    set(H.OnOff_mode,'Value',1);
                    disp('Onset Mode set to 1 - source 1 alone uses random onsets and offsets');
                    [source1_L] = ramp_sound_onoff(source1_L,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                    [source1_R] = ramp_sound_onoff(source1_R,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
            end
            %end
        
        % remove any DCoffset
        source1_L = source1_L - mom(source1_L,1);
        source1_R = source1_R - mom(source1_R,1);
        
        % replace deleted points with zeros
        source1_L_temp=zeros(1,origPnts);
        source1_L_temp(-Delay_pnts:origPnts) = source1_L;
        source1_L = source1_L_temp;
        
        source1_R_temp=zeros(1,origPnts);
        source1_R_temp(-Delay_pnts:origPnts) = source1_R;
        source1_R = source1_R_temp;
        
    else
    
        %if get(H.stim_type,'Value') ~= 9      %Ramp the stimulus
            %ramp_time_on = 10; %ms
            %ramp_time_off = 2.5; %ms
            switch get(H.DelayOnOff_mode,'Value') % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
                case 1
                    [source1_L] = ramp_sound_onoff(source1_L,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                    [source1_R] = ramp_sound_onoff(source1_R,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                case 2
                    [source1_L] = ramp_sound_onoff(source1_L,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                    [source1_R] = ramp_sound_onoff(source1_R,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                case 3
                    [source1_L] = ramp_sound_onoff(source1_L,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
                    [source1_R] = ramp_sound_onoff(source1_R,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
                otherwise
                    set(H.OnOff_mode,'Value',1);
                    disp('Onset Mode set to 1 - source 1 alone uses random onsets and offsets');
                    [source1_L] = ramp_sound_onoff(source1_L,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                    [source1_R] = ramp_sound_onoff(source1_R,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
            end
        %end
        
            % remove any DCoffset
            source1_L = source1_L - mom(source1_L,1);
            source1_R = source1_R - mom(source1_R,1);
            
    end
    
        
    %Apply ITD filtering if conducting ILDAlone Two_source Test
    if(XStimParams.ildalone_flag == 1)
        itdleft = 0; itdright = 0;
        ITD = round(str2num(get(H.ITD,'String')));
        if(ITD < 0)
            itdleft = 0;
            itdright = abs(ITD);
        elseif(ITD > 0)
            itdleft = abs(ITD);
            itdright = 0;
        end
        if(trialnum == 1)
            eval(['load ' FN.ITD_path 'itdfilt' num2str(itdleft)]);
            eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
            eval(['load ' FN.ITD_path 'itdfilt' num2str(itdright)]);
            eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
        end
        source1_L = conv(source1_L,itd_filt_left);
        source1_R = conv(source1_R,itd_filt_right);
    end
    
    % Apply ILD filtering if conducting ITDalone Test
    if(XStimParams.itdalone_flag == 1)
        if(trialnum == 1)
            %%%%%%%%%%%%%%
            if FN.HRTFfiletype(2) == 1
                dir_ILA = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone]));
                ILAind = max(find(dir_ILA(1,:) == XStimParams.el & dir_ILA(2,:) == XStimParams.az));
                if isempty('ILAind')
                    disp('Could not find ILA location in HRTF file');
                    return
                end
                eval(['ila_filt_left = mtlrch(''' FN.ILA_path FN.ildalone ''', ' num2str(ILAind * 2-1) ');']);
                eval(['ila_filt_right = mtlrch(''' FN.ILA_path FN.ildalone ''', ' num2str(ILAind * 2) ');']);
            else
                ILAind = max(find(dir_ILA(1,:) == XStimParams.el & dir_ILA(2,:) == XStimParams.az));
                if isempty('ILAind')
                    disp('Could not find ILA location in HRTF file');
                    return
                end
                ila_filt_left = TF1_ila(ILAind,:);
                ila_filt_right = TF2_ila(ILAind,:);
            end
            %%%%%%%%%%%%
            
        end
        source1_L = conv(source1_L,ila_filt_left);
        source1_R = conv(source1_R,ila_filt_right);
    end
    
    %Add in the leading silent period
    source1_L =  [zeros(1,silence_len) source1_L];
    source1_R = [zeros(1,silence_len) source1_R];
    
    %Add in the trailing silent period
    source1_L =  [source1_L zeros(1,silence_len2)];
    source1_R = [source1_R zeros(1,silence_len2)];
    
    %Add to account for delays that are shorter than maximum
    if(XStimParams.randOnsetPerms)
        source1_L =  [source1_L zeros(1,MaxDelayBufferPnts)];
        source1_R = [source1_R zeros(1,MaxDelayBufferPnts)];
    end
    
    %Apply HRTF filtering
    if(XStimParams.space_flag == 1)
        if FN.HRTFfiletype(1) == 1
            hrtf_left = mtlrch([FN.space_path FN.space],(2*LocIndSeq(1,randseq(trialnum)))-1); % !! locind -> LocIndSeq
            hrtf_right = mtlrch([FN.space_path FN.space],2*LocIndSeq(1,randseq(trialnum)));    % !! locind -> LocIndSeq
        else
            hrtf_left = TF1_space(LocIndSeq(1,randseq(trialnum)),:);                           % !! locind -> LocIndSeq
            hrtf_right = TF2_space(LocIndSeq(1,randseq(trialnum)),:);                          % !! locind -> LocIndSeq
        end
    elseif(XStimParams.ABLalone_flag == 1)
        % get ABLfactor for this location
        if FN.HRTFfiletype(1) == 1
            tempL = mtlrch([FN.space_path FN.space],(2*LocIndSeq(1,randseq(trialnum)))-1);     % !! locind -> LocIndSeq
            tempR = mtlrch([FN.space_path FN.space],2*LocIndSeq(1,randseq(trialnum)));         % !! locind -> LocIndSeq
        else
            tempL = TF1_space(LocIndSeq(1,randseq(trialnum)),:);                               % !! locind -> LocIndSeq
            tempR = TF2_space(LocIndSeq(1,randseq(trialnum)),:);                               % !! locind -> LocIndSeq
        end
        tempL = ERBFilterBankB(tempL, fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
        tempR = ERBFilterBankB(tempR, fcoefs) .* Factormat;
        [ildx ablx] = calclevel_time(tempL,tempR, cF);
        %ITDx = calcitd(tempL,tempR, cF, Fs, ones(n_cF,1));
        ABLfactor = (10 .^(ablx/20)) ./ (10 .^(ablref/20));
        %%%%%%%%%%%%%%%%%% next line
        ABLfactor_long = (interp1(freq(freq_ind),ABLfactor,freq(startpt:endpt)))';
        % apply ABLfactor to reference location
        FT_L = fft(Lref,fftpts);
        FT_L(Xpart1) = FT_L(Xpart1) .* (ABLfactor_long);       	% positive freqs
        FT_L(Xpart2) = FT_L(Xpart2) .* flipud(ABLfactor_long);   % negative freqs
        hrtf_left = real(ifft(FT_L));
        hrtf_left = hrtf_left(1:255);
        FT_R = fft(Rref,fftpts);
        FT_R(Xpart1)= FT_R(Xpart1) .* (ABLfactor_long);       	% positive freqs
        FT_R(Xpart2)= FT_R(Xpart2) .* flipud(ABLfactor_long);   % negative freqs
        hrtf_right = real(ifft(FT_R));
        hrtf_right = hrtf_right(1:255);
        clear FT_L FT_R ABLfactor*         
    elseif(XStimParams.ildalone_flag == 1)
        if FN.HRTFfiletype(2) == 1
            hrtf_left = mtlrch([FN.ILA_path FN.ildalone],(2*LocIndSeq(1,randseq(trialnum)))-1);    % !! locind -> LocIndSeq
            hrtf_right = mtlrch([FN.ILA_path FN.ildalone],2*LocIndSeq(1,randseq(trialnum)));       % !! locind -> LocIndSeq
        else
            hrtf_left = TF1_ila(LocIndSeq(1,randseq(trialnum)),:);                                 % !! locind -> LocIndSeq
            hrtf_right = TF2_ila(LocIndSeq(1,randseq(trialnum)),:);                                % !! locind -> LocIndSeq
        end
    elseif(XStimParams.itdalone_flag == 1)
        if FN.HRTFfiletype(3) == 1
            hrtf_left = mtlrch([FN.ITA_path FN.itdalone],(2*LocIndSeq(1,randseq(trialnum)))-1);    % !! locind -> LocIndSeq
            hrtf_right = mtlrch([FN.ITA_path FN.itdalone],2*LocIndSeq(1,randseq(trialnum)));       % !! locind -> LocIndSeq
        else
            hrtf_left = TF1_ita(LocIndSeq(1,randseq(trialnum)),:);                                 % !! locind -> LocIndSeq
            hrtf_right = TF2_ita(LocIndSeq(1,randseq(trialnum)),:);                                % !! locind -> LocIndSeq
        end
    end
    source1_L = conv(source1_L,hrtf_left);
    source1_R = conv(source1_R,hrtf_right);
    
    % scale
    source1_L = source1_L - round(mean(source1_L));
    source1_R = source1_R - round(mean(source1_R));
    if XStimParams.ildalone_flag | XStimParams.itdalone_flag
        ABAval = 0.5*(mom(source1_L,2) + mom(source1_R,2));
        scalefact = TDT.scalevalue/ABAval;
        source1_L = round(scalefact*source1_L);
        source1_R = round(scalefact*source1_R);
    end
    
   %%%%%%%%% construct uncorrelated stim OR just filter the 2nd stimulus 
   if(uncorrel) % new stim 2
        if get(H.stim_type,'Value') ~= 9
            source2_L = get_stim(XStimParams.freq(2));
        else
            source2_L = get_stim(FN.stim_path2,FN.stim2);
        end
        source2_R = source2_L; 
        
        % modulate stim2
        if ~strcmp(XStimParams.mod_type,'None')
            Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
            source2_R = source2_R .* Envelope(:)';
            source2_L = source2_L .* Envelope(:)';
        end
    end
    
    % onset type (value for SilentDelay) set above for source1
    if(Delay_pnts>0 & SilentDelay) % zero time of delay in second stim.
        origPnts=size(source2_L,2);
        source2_L_temp = source2_L;
        source2_R_temp = source2_R;
        
        % delete delay
        %disp(size(source1_L_temp,2));
        source2_L = source2_L_temp(Delay_pnts:size(source2_L_temp,2));
        source2_R = source2_R_temp(Delay_pnts:size(source2_R_temp,2));
        %disp(size(source1_L,2));
        
       % if get(H.stim_type,'Value') ~= 9      %Ramp the stimulus
            %ramp_time_on = 10; %ms
            %ramp_time_off = 2.5; %ms
            switch get(H.DelayOnOff_mode,'Value') % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
                case 1
                    [source2_L] = ramp_sound_onoff(source2_L,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
                    [source2_R] = ramp_sound_onoff(source2_R,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
                case 2
                    [source2_L] = ramp_sound_onoff(source2_L,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                    [source2_R] = ramp_sound_onoff(source2_R,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                case 3
                    [source2_L] = ramp_sound_onoff(source2_L,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                    [source2_R] = ramp_sound_onoff(source2_R,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                otherwise
                    set(H.OnOff_mode,'Value',1);
                    disp('Onset Mode set to 1 - source 1 alone uses random onsets and offsets');
                    [source2_L] = ramp_sound_onoff(source2_L,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
                    [source2_R] = ramp_sound_onoff(source2_R,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
            end
       % end
        
        % remove any DCoffset
        source2_L = source2_L - mom(source2_L,1);
        source2_R = source2_R - mom(source2_R,1);
        
        % replace deleted points with zeros
        source2_L_temp=zeros(1,origPnts);
        source2_L_temp(Delay_pnts:origPnts) = source2_L;
        source2_L = source2_L_temp;
        
        source2_R_temp=zeros(1,origPnts);
        source2_R_temp(Delay_pnts:origPnts) = source2_R;
        source2_R = source2_R_temp;
           
    else
        %if get(H.stim_type2,'Value') ~= 9  	%Ramp the stimuli
            %ramp_time_on = 10; %ms
            %ramp_time_off = 2.5; %ms
            switch get(H.DelayOnOff_mode,'Value') % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
                case 1
                    [source2_L] = ramp_sound_onoff(source2_L,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
                    [source2_R] = ramp_sound_onoff(source2_R,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
                case 2
                    [source2_L] = ramp_sound_onoff(source2_L,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                    [source2_R] = ramp_sound_onoff(source2_R,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                case 3
                    [source2_L] = ramp_sound_onoff(source2_L,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                    [source2_R] = ramp_sound_onoff(source2_R,TDT.Fs, OnsetSeq(trialnum),OffsetSeq(trialnum));
                otherwise
                    set(H.OnOff_mode,'Value',1);
                    disp('Onset Mode set to 1 - source 1 alone uses random onsets and offsets');
                    [source2_L] = ramp_sound_onoff(source2_L,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
                    [source2_R] = ramp_sound_onoff(source2_R,TDT.Fs, XStimParams.ramp_timeS_on(1), XStimParams.ramp_timeS_off(1));
            end
        %end
        
        % remove any DCoffset
        source2_L = source2_L - mom(source2_L,1);
        source2_R = source2_R - mom(source2_R,1);
    end  
      
          
    %Apply ITD filtering if conducting ILDAlone Two_source Test
    if(XStimParams.ildalone_flag == 1)
        itdleft = 0; itdright = 0;
        ITD = round(str2num(get(H.ITD,'String')));
        if(ITD < 0)
            itdleft = 0;
            itdright = abs(ITD);
        elseif(ITD > 0)
            itdleft = abs(ITD);
            itdright = 0;
        end
        if(trialnum == 1)
            eval(['load ' FN.ITD_path 'itdfilt' num2str(itdleft)]);
            eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
            eval(['load ' FN.ITD_path 'itdfilt' num2str(itdright)]);
            eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
        end
        source2_L = conv(source2_L,itd_filt_left);
        source2_R = conv(source2_R,itd_filt_right);
    end
    % Apply ILD filtering if conducting ITDalone Test
    if(XStimParams.itdalone_flag == 1)
        if(trialnum == 1)
            ILAind = max(find(dir_ILA(1,:) == XStimParams.el & dir_ILA(2,:) == XStimParams.az));
            if isempty('ILAind')
                disp('Could not find ILA location #2 in HRTF file');
                return
            end
            if FN.HRTFfiletype(2) == 1
                eval(['ila_filt_left2 = mtlrch(''' FN.ILA_path FN.ildalone ''', ' num2str(ILAind * 2-1) ');']);
                eval(['ila_filt_right2 = mtlrch(''' FN.ILA_path FN.ildalone ''', ' num2str(ILAind * 2) ');']);
            else
                ila_filt_left2 = TF1_ila(ILAind,:);
                ila_filt_right2 = TF2_ila(ILAind,:);
            end
        end
        source2_L = conv(source2_L,ila_filt_left2);
        source2_R = conv(source2_R,ila_filt_right2);
    end
        
    %Add in the leading silent period
    source2_L =  [zeros(1,silence_len) source2_L];
    source2_R = [zeros(1,silence_len) source2_R];
    
    %Add in the trailing silent period
    source2_L =  [source2_L zeros(1,silence_len2)];
    source2_R = [source2_R zeros(1,silence_len2)];
    
    %Add to account for delays that are shorter than maximum
    if(XStimParams.randOnsetPerms)
        source2_L =  [source2_L zeros(1,MaxDelayBufferPnts)];
        source2_R = [source2_R zeros(1,MaxDelayBufferPnts)];
    end
    
    %Apply HRTF filtering
    if(XStimParams.space_flag == 1)
        if FN.HRTFfiletype(1) == 1
            hrtf_left = mtlrch([FN.space_path FN.space],(2*LocIndSeq(2,randseq(trialnum)))-1); % !! locind -> LocIndSeq
            hrtf_right = mtlrch([FN.space_path FN.space],2*LocIndSeq(2,randseq(trialnum)));    % !! locind -> LocIndSeq
        else
            hrtf_left = TF1_space(LocIndSeq(2,randseq(trialnum)),:);                           % !! locind -> LocIndSeq
            hrtf_right = TF2_space(LocIndSeq(2,randseq(trialnum)),:);                          % !! locind -> LocIndSeq
        end
    elseif(XStimParams.ABLalone_flag == 1)
        % get ABLfactor for this location
        if FN.HRTFfiletype(1) == 1
            tempL = mtlrch([FN.space_path FN.space],(2*LocIndSeq(2,randseq(trialnum)))-1);     % !! locind -> LocIndSeq
            tempR = mtlrch([FN.space_path FN.space],2*LocIndSeq(2,randseq(trialnum)));         % !! locind -> LocIndSeq
        else
            tempL = TF1_space(LocIndSeq(2,randseq(trialnum)),:);                               % !! locind -> LocIndSeq
            tempR = TF2_space(LocIndSeq(2,randseq(trialnum)),:);                               % !! locind -> LocIndSeq
        end
        tempL = ERBFilterBankB(tempL, fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
        tempR = ERBFilterBankB(tempR, fcoefs) .* Factormat;
        [ildx ablx] = calclevel_time(tempL,tempR, cF);
        ABLfactor = (10 .^(ablx/20)) ./ (10 .^(ablref/20));
        %%%%%%%%%%%%% next line
        ABLfactor_long = (interp1(freq(freq_ind),ABLfactor,freq(startpt:endpt)))';
        % apply ABLfactor to reference location
        FT_L = fft(Lref,fftpts);
        FT_L(Xpart1)= FT_L(Xpart1) .* (ABLfactor_long);       	% positive freqs
        FT_L(Xpart2)= FT_L(Xpart2) .* flipud(ABLfactor_long);   % negative freqs
        hrtf_left = real(ifft(FT_L));
        hrtf_left = hrtf_left(1:255);
        FT_R = fft(Rref,fftpts);
        FT_R(Xpart1)= FT_R(Xpart1) .* (ABLfactor_long);       	% positive freqs
        FT_R(Xpart2)= FT_R(Xpart2) .* flipud(ABLfactor_long);   % negative freqs
        hrtf_right = real(ifft(FT_R));
        hrtf_right = hrtf_right(1:255);
        clear FT_L FT_R ABLfactor* 
    elseif(XStimParams.ildalone_flag == 1)
        if FN.HRTFfiletype(2) == 1
            hrtf_left = mtlrch([FN.ILA_path FN.ildalone],(2*LocIndSeq(2,randseq(trialnum)))-1);    % !! locind -> LocIndSeq
            hrtf_right = mtlrch([FN.ILA_path FN.ildalone],2*LocIndSeq(2,randseq(trialnum)));       % !! locind -> LocIndSeq
        else
            hrtf_left = TF1_ila(LocIndSeq(2,randseq(trialnum)),:);                                 % !! locind -> LocIndSeq
            hrtf_right = TF2_ila(LocIndSeq(2,randseq(trialnum)),:);                                % !! locind -> LocIndSeq
        end
    elseif(XStimParams.itdalone_flag == 1)
        if FN.HRTFfiletype(3) == 1
            hrtf_left = mtlrch([FN.ITA_path FN.itdalone],(2*LocIndSeq(2,randseq(trialnum)))-1);    % !! locind -> LocIndSeq
            hrtf_right = mtlrch([FN.ITA_path FN.itdalone],2*LocIndSeq(2,randseq(trialnum)));       % !! locind -> LocIndSeq
        else
            hrtf_left = TF1_ita(LocIndSeq(2,randseq(trialnum)),:);                             % !! locind -> LocIndSeq
            hrtf_right = TF2_ita(LocIndSeq(2,randseq(trialnum)),:);                            % !! locind -> LocIndSeq
        end
    end
    source2_L = conv(source2_L,hrtf_left);
    source2_R = conv(source2_R,hrtf_right);
    
    % scale
    source2_L = source2_L - round(mean(source2_L));
    source2_R = source2_R - round(mean(source2_R));
    if XStimParams.ildalone_flag | XStimParams.itdalone_flag
        ABAval = 0.5*(mom(source2_L,2) + mom(source2_R,2));
        scalefact = TDT.scalevalue/ABAval;
        source2_L = round(scalefact*source2_L);
        source2_R = round(scalefact*source2_R);
    end
    
    
    %source2_L = source1_L;
    %source2_R = source1_R;
    
    % Delay
    if(Delay_pnts>0) % delay 1st sound
        source1_L_temp=zeros(1,size(source1_L,2)+Delay_pnts);
        source1_L_temp(size(source1_L_temp,2)-size(source1_L,2)+1:size(source1_L_temp,2))=source1_L;
        source1_L = source1_L_temp;
        
        source1_R_temp=zeros(1,size(source1_R,2)+Delay_pnts);
        source1_R_temp(size(source1_R_temp,2)-size(source1_R,2)+1:size(source1_R_temp,2))=source1_R;
        source1_R = source1_R_temp;
        
        
        source2_R_temp=zeros(1,size(source2_R,2)+Delay_pnts);
        source2_R_temp(1:size(source2_R,2))=source2_R;
        source2_R = source2_R_temp;
        
        source2_L_temp=zeros(1,size(source2_L,2)+Delay_pnts);
        source2_L_temp(1:size(source2_L,2))=source2_L;
        source2_L = source2_L_temp;
        
    else
        if(Delay_pnts<0) % delay 2nd sound
            source2_L_temp=zeros(1,size(source2_L,2)-Delay_pnts);
            source2_L_temp(size(source2_L_temp,2)-size(source2_L,2)+1:size(source2_L_temp,2))=source2_L;
            source2_L = source2_L_temp;
            
            source2_R_temp=zeros(1,size(source2_R,2)-Delay_pnts);
            source2_R_temp(size(source2_R_temp,2)-size(source2_R,2)+1:size(source2_R_temp,2))=source2_R;
            source2_R = source2_R_temp;
            
            
            source1_R_temp=zeros(1,size(source2_R,2)-Delay_pnts);
            source1_R_temp(1:size(source2_R,2))=source2_R;
            source1_R = source2_R_temp;
            
            source1_L_temp=zeros(1,size(source1_L,2)-Delay_pnts);
            source1_L_temp(1:size(source1_L,2))=source1_L;
            source1_L = source1_L_temp;
        end
    end
    
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
    
    %if(1 |trialnum==1)
    %if(DelayModsSeq(trialnum)==4)
    if(1 |trialnum==1)
        timeaxis = 0 : 1/TDT.Fs : (size(source1_L,2)*(1/TDT.Fs))-(1/TDT.Fs);
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
        Stim.L = plot(timeaxis, source1_L + 0.4,...
            'Tag','Env_ellipses',...
            'Marker', 'none',...
            'LineStyle','-',...
            'Color', [1 0 0],...
            'LineWidth', 0.5);
        %subplot(2,1,2);
        Stim.R = plot(timeaxis, source2_L - 0.4,...
            'Tag','Env_ellipses',...
            'Marker', 'none',...
            'LineStyle','-',...
            'Color', [0 0 1],...
            'LineWidth', 0.5);
        hold off;
        %set(gca, 'Title', 'red=src1_L, blue=src2_L');
        title('red=src1_L, blue=src2_L');
    end

    % save stims to disk with name of loc1
    % LocationsSeq, LocIndSeq, DelayTimesSeq, OnsetSeq, OffsetSeq, LagExtSegsSeq, DelayModsSeq
    if(exist1('H.Delayfig'));
        S232('push16',filttrial_left,length(filttrial_left));
        S232('qpop16',BUF.L1);
        fname = ['stimbuf_left_' num2str(hrtfdirmat(1,LocIndSeq(1,randseq(trialnum)))) ...
                '_' num2str(hrtfdirmat(2,LocIndSeq(1,randseq(trialnum)))) ...
                '_' num2str(DelayTimesSeq(randseq(trialnum))) ...
                '_' num2str(OnsetSeq(randseq(trialnum))) ...
                '_' num2str(OffsetSeq(randseq(trialnum))) ...
                '_' num2str(LagExtSegsSeq(randseq(trialnum))) ...
                '_' num2str(DelayModsSeq(randseq(trialnum)))];        % !! 2X locind -> LocIndSeq
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
                '_' num2str(DelayModsSeq(randseq(trialnum)))];        % !! 2X locind -> LocIndSeq
        %disp(fname);
        evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
                [' ''' stimuli_dir fname ''' ']   ',0);'];
        eval(evalstr);
        temp_right = dama2pc(BUF.R1);
    end
    
    remtrials = numtrials - trialnum;
    set(H.remtrials,'String',num2str(remtrials));
    trialnum = trialnum + 1;
    set(H.status,'BackgroundColor','blue');
    pause(0);
end 										%end loop over trials

%%%%%%%%%%%%%%%%%%%%%%%%%%% finished making sounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin playing sounds   
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
    %Randomize the stimuli
    randseq = randperm(numtrials);   
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
    
    % loop for trials
    tic
    while (exist1('H.Delayfig') & (trialnum <= numtrials+1))
        
        %Check for pause by user
        if pause_check  return; end
        
        %Wait till PD1 is finished
        while S232('PD1status',1) usec_delay(1000); end
        
        % LocationsSeq, LocIndSeq, DelayTimesSeq, OnsetSeq, OffsetSeq, LagExtSegsSeq, DelayModsSeq
        if(trialnum <= numtrials)
            fname = ['stimbuf_left_' num2str(hrtfdirmat(1,LocIndSeq(1,randseq(trialnum)))) ...
                    '_' num2str(hrtfdirmat(2,LocIndSeq(1,randseq(trialnum)))) ...
                    '_' num2str(DelayTimesSeq(randseq(trialnum))) ...
                    '_' num2str(OnsetSeq(randseq(trialnum))) ...
                    '_' num2str(OffsetSeq(randseq(trialnum))) ...
                    '_' num2str(LagExtSegsSeq(randseq(trialnum))) ...
                    '_' num2str(DelayModsSeq(randseq(trialnum)))];                                   % !! 2X locind -> LocIndSeq%
            evalstr = ['S232(''disk2dama16'',BUF.L1,'  [' ''' stimuli_dir fname ''' '] ',0);'];
            eval(evalstr);
            fname = ['stimbuf_right_' num2str(hrtfdirmat(1,LocIndSeq(1,randseq(trialnum)))) ...
                    '_' num2str(hrtfdirmat(2,LocIndSeq(1,randseq(trialnum)))) ...
                    '_' num2str(DelayTimesSeq(randseq(trialnum))) ...
                    '_' num2str(OnsetSeq(randseq(trialnum))) ...
                    '_' num2str(OffsetSeq(randseq(trialnum))) ...
                    '_' num2str(LagExtSegsSeq(randseq(trialnum))) ...
                    '_' num2str(DelayModsSeq(randseq(trialnum)))];                                   % !! 2X locind -> LocIndSeq%
            evalstr = ['S232(''disk2dama16'',BUF.R1,' [' ''' stimuli_dir fname ''' '] ',0);'];
            eval(evalstr);
        end
        
        %Wait till PD1 is finished
        while S232('PD1status',1) usec_delay(1000); end
        S232('PD1stop',1);
        
        %Stop the m110 and get spikes
        if(trialnum > 1)							% first trial just for loading
            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, 2*XStimParams.curr_stimdur); 			% Take XStimParams.curr_stimdur spikes max
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
              'Position',[700 100 400 400],...
              'Tag', 'Raster',...
              'Name','Raster',...
              'NumberTitle','off',...
              'Color',[1 1 1]);
            H.delay_raster_axes = axes;
            figure(H.delay_raster);
            ylabel('trial'); xlabel('time (ms)'); title(['Rep# ' num2str(repnum)]);
            ylim([0 numtrials+1]);
            set(H.delay_raster_axes, 'YLim', [0 numtrials+0.5]);
        end
        plotmatrix = [Nspikes spikes_trial repnum_trial EL_trial AZ_trial Delay_trial  Onset_ms_trial Offset_ms_trial LagExtSegs_ms_trial OnsetMod_trial];
        %assignin('base','plotmatrix',plotmatrix);
        Delay_raster(XStimParams,plotmatrix, randseq, XStimParams.numreps);
        
        
        %Record Data
        if(exist1('H.Delayfig') & get(H.recorddata,'Value') )
            tempseq{repnum} = randseq;
            datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial Delay_trial  Onset_ms_trial Offset_ms_trial LagExtSegs_ms_trial OnsetMod_trial]];
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

% % plot results over all reps
% if(exist1('H.Delayfig') & get(H.recorddata,'Value') )
% 	if(findobj('Tag','Raster'))
%         close('Raster');
% 	end
% 	H.delay_raster = figure('Units','pixels',...
%       'Position',[700 100 400 400],...
%       'Tag', 'Raster',...
%       'Name','Raster',...
%       'NumberTitle','off',...
%       'Color',[1 1 1]);
% 	H.delay_raster_axes = axes;
% 	figure(H.delay_raster);
% 	ylabel('trial'); xlabel('time (ms)'); title([num2str(repnum-1) ' Reps ']);
% 	ylim([0 numtrials+1]);
% 	set(H.delay_raster_axes, 'YLim', [0 numtrials]);
% 	
% 	Delay_raster(XStimParams,datamatrix, randseq);
% end

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

%/////////////////////////////////////////////////////////////
%%%%%%%%%
function [stim] = get_stim(param1, param2)
% param1: XStimParams.freq or FN.stim_path
% param2: FN.stim

global H
global XStimParams
global TDT

switch get(H.stim_type,'Value')      
    case 1    %tone at specified frequency 
        stim = MakeTone(TDT.Fs,param1,XStimParams.curr_stimdur+50);
    case 2    %GammaTones
        stim = MakeGammaTone(TDT.Fs,param1,XStimParams.curr_stimdur+50);
    case 8    %Broadband Noise
        [stim] = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur+50);
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
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur+50);
        disp('Stimulus type not supported for Delay tests.  Reset to BROADBAND');
        return
end

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