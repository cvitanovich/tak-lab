function [] = Engage_abltest()

global H
global XStimParams
global TDT
global FN
global C_
global M_

%Engage_abltest

%*******************************************************************************
%	The abltest operation
%*******************************************************************************
stimuli_dir = FN.temp_stim_path;
fclose all;
if ~XStimParams.buildOnly
    eval(['delete ' stimuli_dir '*.*;']);
end
disp('This is an ABL test')

% force use of either BBnoise or tone
if get(H.stim_type,'Value') ~= 1 & get(H.stim_type,'Value') ~= 8;
    set(H.stim_type,'Value',8);         % broadband noise
end

% use HRTFs etc.
UseFilts = 0;                           % itd only
if XStimParams.ephone_flag
    UseFilts = 1;                   % itd and ephone
end
UseFilts = XStimParams.loc_flag + UseFilts;     % std or eq hrtfs

%Put parameters into XStimParams
XStimParams.test_type = 'ABL';
XStimParams.loabl = str2num(get(H.lowabl,'String'));
XStimParams.hiabl = str2num(get(H.highabl,'String'));
XStimParams.numabls = str2num(get(H.numabls,'String'));
XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
XStimParams.curr_ILD = str2num(get(H.ILD,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.reset_flag = 0;

Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_ABL Temp_params;'])
clear Temp_params

%Setup for ABLs used
abls  = XStimParams.abls;

numreps = XStimParams.numreps;
numtrials = length(abls);

clear BUF

%Specify DAMA buffers
BUF.left_1			= 1;
BUF.right_1		= 2;
BUF.playseq_left_1			= 4;
BUF.playseq_right_1		= 5;
BUF.playspecbuf_1			= 6;
BUF.left_2			= 7;
BUF.right_2		= 8;
BUF.playseq_left_2			= 9;
BUF.playseq_right_2		= 10;
BUF.playspecbuf_2			= 11;

%Make play sequence buffers
S232('allot16',BUF.playseq_left_1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.left_1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_left_1);

S232('allot16',BUF.playseq_right_1,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.right_1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_right_1);

S232('allot16',BUF.playseq_left_2,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.left_2);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_left_2);

S232('allot16',BUF.playseq_right_2,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.right_2);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_right_2);

%Make play specification buffer
S232('allot16',BUF.playspecbuf_1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_left_1);
S232('make',1,BUF.playseq_right_1);
S232('make',2,0);
S232('qpop16',BUF.playspecbuf_1);

S232('allot16',BUF.playspecbuf_2,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_left_2);
S232('make',1,BUF.playseq_right_2);
S232('make',2,0);
S232('qpop16',BUF.playspecbuf_2);

%load ephone or ephone2 files
str2 = 'no hrtfs';
switch UseFilts
    case 1              % ephone
        if isempty(FN.ephone)        % FN not yet picked
            [FN.ephone,FN.ephone_path] = uigetfile('*.*','Select ephone files for use withOUT HRTFs (onesX.imp preferred)');
            FN.HRTFfiletype(5,1) = testHRTFfiletype(FN.ephone_path, FN.ephone);
        end
        if FN.HRTFfiletype(5,1) == 1
            ephonefilt_left = mtlrch([FN.ephone_path FN.ephone],1);
            ephonefilt_right = mtlrch([FN.ephone_path FN.ephone],2);
        elseif FN.HRTFfiletype(5,1) == 2
            dir = 0;
            eval(['load -mat ' FN.ephone_path FN.ephone]);
            ephonefilt_left = TF1(1,:);
            ephonefilt_right = TF2(1,:);
            clear dir TF1 TF2
        else
            disp(['ephone files incorrect'])
            return
        end
    case 2              % ephone2 and std hrtfs 
        if isempty(FN.ephone2)       % FN not yet picked
            [FN.ephone2,FN.ephone_path] = uigetfile('*.*','Select ephone2 files for use with HRTFs (ET_EC_inv preferred)');
            FN.HRTFfiletype(6,1) = testHRTFfiletype(FN.ephone_path, FN.ephone2);
        end
        if FN.HRTFfiletype(6,1) == 1
            ephonefilt_left = mtlrch([FN.ephone_path FN.ephone2],1);
            ephonefilt_right = mtlrch([FN.ephone_path FN.ephone2],2);
        elseif FN.HRTFfiletype(6,1) == 2
            dir = 0;
            eval(['load -mat ' FN.ephone_path FN.ephone2]);
            ephonefilt_left = TF1(1,:);
            ephonefilt_right = TF2(1,:);
            clear dir TF1 TF2
        else
            disp(['ephone2 files incorrect'])
            return
        end
end

% load HRTFs
switch UseFilts
    case 2
        if isempty(FN.space_std)       % FN not yet picked
            [FN.space_std,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.std preferred)');
            FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
        end
        if FN.HRTFfiletype(1,2) == 1
            dirmat = sph2dbl(mtlrdir([FN.space_path FN.space_std]));
        elseif FN.HRTFfiletype(1,2) == 2
            dir = 0;
            eval(['load -mat ' FN.space_path FN.space_std]);
            dirmat = dir;
            clear dir
        else
            disp(['space_std HRTFfiletype incorrect'])
            return
        end
        locind = max(find(dirmat(1,:) == XStimParams.loc_azel(2) & dirmat(2,:) == XStimParams.loc_azel(1)));
        % get HRTFs
        if FN.HRTFfiletype(1,2) == 1
            hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locind)-1);
            hrtf_right = mtlrch([FN.space_path FN.space_std],2*locind);
        else
            hrtf_left = TF1(locind,:);
            hrtf_right = TF2(locind,:);
            clear TF1 TF2
        end
        str2 = ['with hrtf.std at ' num2str(XStimParams.loc_azel(2)) ' el, ' num2str(XStimParams.loc_azel(1)) ' az'];
    case 3
        if isempty(FN.space_eq)       % FN not yet picked
            [FN.space_eq,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.eq preferred)');
            FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);
        end
        if FN.HRTFfiletype(1,1) == 1
            dirmat = sph2dbl(mtlrdir([FN.space_path FN.space_eq]));
        elseif FN.HRTFfiletype(1,1) == 2
            dir = 0;
            eval(['load -mat ' FN.space_path FN.space_eq]);
            dirmat = dir;
            clear dir
        else
            disp(['space_eq HRTFfiletype incorrect'])
            return
        end
        locind = max(find(dirmat(1,:) == XStimParams.loc_azel(2) & dirmat(2,:) == XStimParams.loc_azel(1)));
        % get HRTFs
        if FN.HRTFfiletype(1,1) == 1
            hrtf_left = mtlrch([FN.space_path FN.space_eq],(2*locind)-1);
            hrtf_right = mtlrch([FN.space_path FN.space_eq],2*locind);
        else
            hrtf_left = TF1(locind,:);
            hrtf_right = TF2(locind,:);
            clear TF1 TF2
        end
        str2 = ['with hrtf.eq at ' num2str(XStimParams.loc_azel(2)) ' el, ' num2str(XStimParams.loc_azel(1)) ' az'];
end

% write diary
if get(H.recorddata,'Value')
    tempstr = ['    ABL-test: ' str2 ];
    update_diary
end

% Add a piece of silence prior to stimulus to calculate spontaneous rate
DUR_silence = XStimParams.silence_lead; 					%ms
silence_len = (DUR_silence * round(TDT.Fs/1000));
% Add a piece of silence after stimulus 
DUR_silence2 = XStimParams.silence_trail; 					%ms
silence_len2 = (DUR_silence2 * round(TDT.Fs/1000));

%Make Stimulus buffers
DUR = str2num(get(H.DUR,'String'));
switch UseFilts
    case 0
        Npts_total_play = silence_len + DUR*round(TDT.Fs/1000) + TDT.itdfiltlen;
    case 1
        Npts_total_play = silence_len + DUR*round(TDT.Fs/1000) + TDT.itdfiltlen + TDT.ephonefiltlen;
    case 2
        Npts_total_play = silence_len + DUR*round(TDT.Fs/1000) + TDT.ephonefiltlen + TDT.hrtffiltlen;
    case 3
        Npts_total_play = silence_len + DUR*round(TDT.Fs/1000) + TDT.hrtffiltlen;       
end

S232('allot16',BUF.left_1,Npts_total_play);
S232('allot16',BUF.right_1,Npts_total_play);
S232('allot16',BUF.left_2,Npts_total_play);
S232('allot16',BUF.right_2,Npts_total_play);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,Npts_total_play);

%Load Earphone filters
switch UseFilts
    case 0          % no ephones, no HRTFs    
        S232('PD1clrsched',1);
        S232('PD1nstrms',1,2,0);
        S232('PD1specIB',1,S232('IB',0),S232('DAC',0));
        S232('PD1specIB',1,S232('IB',1),S232('DAC',1));
        
    case 1      %Use ephones
        dspid_left = 0; dspid_right = 1;
        S232('PD1clrsched',1);
        S232('PD1nstrms',1,2,0);
        S232('PD1resetDSP',1,hex2dec('FFF'));
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
        
    case 2          % Use ephone2 and HRTFs
        dspid_left_loc = 0; dspid_right_loc = 1;
        dspid_left_ephone = 2; dspid_right_ephone = 3;
        S232('PD1clrsched',1);
        S232('PD1nstrms',1,2,0);
        S232('PD1resetDSP',1,hex2dec('FFF'));
        S232('dropall');
        %Make connections for left ear
        S232('PD1addsimp',1,S232('DSPout',dspid_left_ephone),S232('DAC',0)); %ephone to DAC0
        S232('PD1addsimp',1,S232('DSPout',dspid_left_loc),S232('DSPin',dspid_left_ephone)); %loc to ephone
        S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left_loc)); %IB to loc
        %Make connections for right ear
        S232('PD1addsimp',1,S232('DSPout',dspid_right_ephone),S232('DAC',1)); %ephone to DAC1
        S232('PD1addsimp',1,S232('DSPout',dspid_right_loc),S232('DSPin',dspid_right_ephone)); %loc to ephone
        S232('PD1specIB',1,S232('IB',1),S232('DSPin',dspid_right_loc)); %IB to loc
        
        %Load left      
        S232('pushf',ephonefilt_left,length(ephonefilt_left));
        S232('PreLoadRaw',1,S232('DSPid',dspid_left_ephone),'MONO','STACK','','',TDT.ephonescale,1.0,1);
        S232('dropall');
        S232('pushf',hrtf_left,length(hrtf_left));
        S232('PreLoadRaw',1,S232('DSPid',dspid_left_loc),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
        S232('dropall');
        %Load right
        S232('pushf',ephonefilt_right,length(ephonefilt_right));
        S232('PreLoadRaw',1,S232('DSPid',dspid_right_ephone),'MONO','STACK','','',TDT.ephonescale,1.0,1);
        S232('dropall');
        S232('pushf',hrtf_right,length(hrtf_right));
        S232('PreLoadRaw',1,S232('DSPid',dspid_right_loc),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
        S232('dropall');
        
    case 3      %Use HRTFs only
        dspid_left = 0; dspid_right = 1;
        S232('PD1clrsched',1);
        S232('PD1nstrms',1,2,0);
        S232('PD1resetDSP',1,hex2dec('FFF'));
        S232('dropall');
        %Make connections for left ear
        S232('PD1addsimp',1,S232('DSPout',dspid_left),S232('DAC',0)); %DSPout to DAC0
        S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left)); %IB to DSPin
        %Make connections for right ear
        S232('PD1addsimp',1,S232('DSPout',dspid_right),S232('DAC',1));
        S232('PD1specIB',1,S232('IB',1),S232('DSPin',dspid_right));
        %Load left      
        S232('pushf',hrtf_left,length(hrtf_left));
        S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
        %Load right
        S232('pushf',hrtf_right,length(hrtf_right));
        S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
end

set(H.ephonefile,'Enable','off');
set(H.ephoneuseit,'Enable','off');
set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

%Set MII parameters
mii_us_per_sample = 10; %microsecond per sample
mii_separation = 100; %only take events separated by 100 samples (i.e., 1 ms)

ITD = XStimParams.curr_ITD;
if(abs(ITD) > 250) return; end

ISI = XStimParams.test_ISI;
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000)); %correct for ITD filtlength

% make the stimuli we'll use
set(H.buildplay,'String','Building Stimuli');
set(H.remreps,'String',num2str(1));
buffcycle = 1;
numtrials = length(abls);
%Randomize the stimuli
randseq = randperm(numtrials);

%Loop to make stimuli
trialnum = 1;
while (exist1('H.ablfig') & (trialnum <= numtrials))
    set(H.buildplay,'BackgroundColor','red');
    %Check for pause by user
    if pause_check  return; end
    
    set(H.stim_type,'Value',8);
    
    %Make at specified frequency and ild - NOT GETTING ILD by ATTENUATORS
    ILD = XStimParams.curr_ILD;
    stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
    
    % remove any DCoffset
    stim = stim - mom(stim,1);
    
    % normalize to ACPower (added 3/5/07 as in 2-source)
    stim = stim / mom(stim,2);
    
    % apply ILD half to each side
    if UseFilts < 2
        trial_left  = stim * (10^(-ILD/(2*20)));
        trial_right = stim * (10^(ILD/(2*20)));
    else
        trial_left = stim;
        trial_right = stim;
    end

    % 3/7/07 multiply by freefield_Xfactors
    if UseFilts < 2
        trial_left = trial_left * TDT.freefield_Lfactor;
        trial_right = trial_right * TDT.freefield_Rfactor;
    end
    
    %Ramp the stimuli
    ramp_time = 5; %ms
    [trial_left] = ramp_sound(trial_left,TDT.Fs,ramp_time);
    [trial_right] = ramp_sound(trial_right,TDT.Fs,ramp_time);
    
    % remove DCoffset
    trial_left = trial_left - mom(trial_left,1);
    trial_right = trial_right - mom(trial_right,1);
    
    %Add in the leading silent period
    trial_left =  [zeros(1,silence_len) trial_left];
    trial_right = [zeros(1,silence_len) trial_right];
    
    %Add in the trailing silent period
    trial_left =  [trial_left zeros(1,silence_len2)];
    trial_right = [trial_right zeros(1,silence_len2)];
    
    %Apply ITD filtering
    if UseFilts < 2
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
            eval(['load ' FN.itd_path 'itdfilt' num2str(itdleft)]);
            eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
            itd_filt_left = itd_filt_left/max(abs(itd_filt_left));
            eval(['load ' FN.itd_path 'itdfilt' num2str(itdright)]);
            eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
            itd_filt_right = itd_filt_right/max(abs(itd_filt_left));
        end
        trial_left = conv(trial_left,itd_filt_left);
        trial_right = conv(trial_right,itd_filt_right);
        % remove DCoffset
        trial_left = trial_left - round(mean(trial_left));
        trial_right = trial_right - round(mean(trial_right));
    end
    
    switch UseFilts
        case 2
            trial_left = [trial_left zeros(1,TDT.ephonefiltlen) zeros(1,TDT.hrtffiltlen)];
            trial_right = [trial_right zeros(1,TDT.ephonefiltlen) zeros(1,TDT.hrtffiltlen)];
        case 3
            trial_left = [trial_left zeros(1,TDT.hrtffiltlen)];
            trial_right = [trial_right zeros(1,TDT.hrtffiltlen)];
    end
    
    % 3/7/07
    trial_left = trial_left * TDT.scaleFactor;
    trial_right = trial_right * TDT.scaleFactor;
    
    % write files
    S232('push16',trial_left,length(trial_left));
    S232('qpop16',BUF.left_1);
    fname = ['ABL.left_' num2str(abls(randseq(trialnum)))];
    evalstr = ['S232(''dama2disk16'',BUF.left_1,' ...
            [' ''' FN.temp_stim_path fname ''' '] ',0);'];
    eval(evalstr);
    temp_left = dama2pc(BUF.left_1);
    S232('push16',trial_right,length(trial_right));
    S232('qpop16',BUF.right_1);
    fname = ['ABL.right_' num2str(abls(randseq(trialnum)))];
    evalstr = ['S232(''dama2disk16'',BUF.right_1,' ...
            [' ''' FN.temp_stim_path fname ''' '] ',0);'];
    eval(evalstr);
    temp_right = dama2pc(BUF.right_1);
    
    %Plot PSD
    if(exist1('H.ablfig') & ~exist1('H.psdaxes'))
        figure(H.ablfig);
        H.psdaxes = axes;
        set(H.psdaxes,'Visible','off');
    end
    if(exist1('H.ablfig') & get(H.plotpsd,'Value') == 1)
        axes(H.psdaxes);
        set(H.psdaxes,'Visible','on',...
            'Position',[0.1 0.07 0.8 0.2]);
        axis square
        [pl,f] = psd(temp_left,2048,TDT.Fs); [pr,f] = psd(temp_right,2048,TDT.Fs);
        hpl = plot(f,10*log10(abs(pl)),'g');
        hold on
        H.pr = plot(f,10*log10(abs(pr)),'r');
        xlabel('Frequency (Hz)');
        ILDval = (20*log10(max(abs(temp_right)))) - (20*log10(max(abs(temp_left))));
        ABLval = 0.5*((20*log10(max(abs(temp_left)))) + (20*log10(max(abs(temp_right)))));
        title(['ABL = ' num2str(ABLval) ', ILD = ' num2str(ILDval)]);
        set(H.psdaxes,'YLim',[0 130]);
        hold off
        grid on
        pause(1)
    elseif(exist1('H.ablfig') & exist1('H.pl') & get(H.plotpsd,'Value') == 0)
        set(H.psdaxes,'Visible','off');
        set(H.pl,'Visible','off');
        set(H.pr,'Visible','off');
    end
    set(H.remtrials,'String',num2str(numtrials-trialnum));
    trialnum = trialnum + 1;
    pause(0);
    
end %end loop over trials

if ~XStimParams.buildOnly
    %Put up raster plot
    position = [1000 0 200 950];
    xvals = 1;
    yvals = abls;
    xtext = '';
    ytext = 'ABL (dB)';
    xoffset = 1;
    yoffset = 1;
    tot_dur = DUR_silence + DUR + 50;
    stim_dur = DUR_silence + DUR;
    startcount = DUR_silence;
    endcount = DUR_silence + DUR;
    row_space = 5;
    col_space = .2*tot_dur;
    
    [H.rasterfig,H.rastaxes] = makerasterfig(position,...
        xvals,...
        yvals,...
        xtext,...
        ytext,...
        xoffset,...
        yoffset,...
        numreps,...
        tot_dur,...
        stim_dur,...
        startcount,...
        endcount,...
        row_space,...
        col_space);
    
    %%%% Begin loop over experimental blocks (i.e., reps) and trials   
    set(H.buildplay,'String','Playing Stimuli');
    set(H.buildplay,'BackgroundColor','yellow');
    set(H.remreps,'String',num2str(numreps));
    repnum = 1;
    buffcycle = 1;
    datamatrix = [];
    finalspikematrix = zeros(1,length(abls));
    
    if(exist1('H.ablfig') & get(H.recorddata,'Value'))
        update_dataFN;
    end
    while (exist1('H.ablfig') & (repnum <= numreps))
        %Randomize the stimuli
        randseq = randperm(numtrials);   
        trialnum = 1;
        spikes_trial = [];
        abls_trial = [];
        repnum_trial = [];
        Nspikes = [];
        
        while (exist1('H.ablfig') & (trialnum <= numtrials+1))
            
            if(trialnum <= numtrials)
                if(buffcycle == 1);
                    fname = ['ABL.left_' num2str(abls(randseq(trialnum)))];
                    evalstr = ['S232(''disk2dama16'',BUF.left_1,' ...
                            [' ''' FN.temp_stim_path fname ''' ']  ',0);'];
                    eval(evalstr);
                    fname = ['ABL.right_' num2str(abls(randseq(trialnum)))];
                    evalstr = ['S232(''disk2dama16'',BUF.right_1,' ...
                            [' ''' FN.temp_stim_path fname ''' ']  ',0);'];
                    eval(evalstr);
                elseif(buffcycle == 2);
                    fname = ['ABL.left_' num2str(abls(randseq(trialnum)))];
                    evalstr = ['S232(''disk2dama16'',BUF.left_2,' ...
                            [' ''' FN.temp_stim_path fname ''' ']  ',0);'];
                    eval(evalstr);
                    fname = ['ABL.right_' num2str(abls(randseq(trialnum)))];
                    evalstr = ['S232(''disk2dama16'',BUF.right_2,' ...
                            [' ''' FN.temp_stim_path fname ''' '] ',0);'];
                    eval(evalstr);
                end
            end
            
            %Wait till PD1 is finished
            %while S232('PD1status',1) usec_delay(1000); end
            usec_delay((DUR+100)*1000)
            S232('PD1stop',1);
            %Check for pause by user
            if pause_check  return;   end
            
            %Stop the m110 and get spikes
            if(trialnum > 1)
                m110dx( C_.STOP);
                spikes = m110dx( C_.DATA, 1000); %Take 1000 spikes max
                ind = find(spikes ~= 0); %Get clock events that are spikes
                spikes = spikes(ind);
                ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
                if(exist1('H.ablfig') & ~isempty(spikes)) 
                    spikes = spikes(ind);
                    spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
                    abls_trial = [abls_trial;randseq(trialnum-1) * ones(size(spikes))];
                    repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
                    Nspikes = [Nspikes;length(spikes) * ones(size(spikes))];
                end
            end
            
            if(exist1('H.ablfig') & trialnum <= numtrials)
                
                if(exist1('H.ablfig') & buffcycle == 1)
                    S232('seqplay',BUF.playspecbuf_1);
                elseif(exist1('H.ablfig') & buffcycle == 2)
                    S232('seqplay',BUF.playspecbuf_2);
                end
                S232('PA4atten',1,abs(abls(randseq(trialnum)))-20);
                S232('PA4atten',2,abs(abls(randseq(trialnum)))-20);
                S232('PD1arm',1);
                
                %Send trigger
                %Set up MII
                m100x( C_.INIT );
                m110dx( C_.INIT );
                m110dx( C_.CLOCK, mii_us_per_sample);
                m110dx( C_.MODE, M_.PST );
                
                pause(ISI/1000);
                if (trialnum <= numtrials)
                    %Start clock
                    m110dx( C_.START);
                    %Send pulse: PD1 GO!
                    m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
                end
            end
            
            if(buffcycle == 1)
                buffcycle = 2;
            else
                buffcycle = 1;
            end
            
            if(trialnum > 1)
                finalspikematrix(randseq(trialnum-1)) = ...
                    finalspikematrix(randseq(trialnum-1)) + ...
                    length(spikes);
            end
            
            remtrials = numtrials - trialnum;
            set(H.remtrials,'String',num2str(remtrials));
            trialnum = trialnum + 1;
            pause(0);
        end %end loop over trials
        
        %Plot all the spikes from trials
        plotraster(H.rasterfig,...
            spikes_trial,...
            ones(size(abls_trial)),...
            abls_trial,...
            repnum,...
            tot_dur,...
            stim_dur,...
            numreps,...
            xoffset,...
            yoffset,...
            row_space,...
            col_space);
        
        %Record Data
        if(exist1('H.ablfig') & get(H.recorddata,'Value'))
            datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial abls_trial]];
            record_data3(XStimParams,datamatrix);
        end
        
        remreps = numreps - repnum;
        set(H.remreps,'String',num2str(remreps));
        repnum = repnum + 1;
        pause(0);
    end %end loop over reps
    
    %Plot final spike rate figure
    finalspikematrix = finalspikematrix/numreps;
    
    H.finalspikeratefig = figure('Position',[700 20 550 500],...
        'Name','Final Test Results',...
        'NumberTitle','off');
    H.axes = axes;
    plot(abls,finalspikematrix,'g-*','LineWidth',2);
    set(H.axes,'Color','black');
    xlabel('ABL (dB)'); ylabel('Spike Rate (spikes/stim)');
    title('ABL Test');
    
    set(H.buildplay,'String','Build/Play status');
end        % if ~XStimParams.buildOnly

% increment testnum
if(exist1('H.ablfig') & get(H.recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
end

%%%%%%%%%
function [flag] = pause_check

global H
global XStimParams
global C_
global M_

flag = 0;
%Check for pause by user
while (exist1('H.ablfig') & get(H.pauseabltest,'Value'))
    pause(0);
    if(~exist1('H.ablfig')) return; end         
    if(exist1('H.ablfig') & get(H.resetabltest,'Value') == 1)
        set(H.resetabltest,'Value',0);
        set(H.pauseabltest,'Value',0);
        Reset_ABLtest;
        flag = 1;
        return;
    end
end
if XStimParams.reset_flag
    flag = 1;
    XStimParams.reset_flag = 0;
end