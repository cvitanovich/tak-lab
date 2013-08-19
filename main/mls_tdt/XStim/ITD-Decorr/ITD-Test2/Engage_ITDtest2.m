function [] = Engage_itdtest2()

global H
global XStimParams
global TDT
global FN
global C_
global M_

%Engage_itdtest2

%*******************************************************************************
%	The itdtest2 operation
% this version allows different sounds for each ear
%*******************************************************************************
stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);
disp('This is an ITD2 test')

% force use of either BBnoise or tone
if get(H.stim_type,'Value') ~= 9
    disp('THIS VERSION OF ITD TEST OPTIMIZED FOR FILES');
end

%Put parameters into XStimParams
test_val = get(H.test_type,'Value');
test_type = get(H.test_type,'String');
XStimParams.test_type = deblank(test_type(test_val,:));
clear test_type
XStimParams.loitd = str2num(get(H.lowitd,'String'));
XStimParams.hiitd = str2num(get(H.highitd,'String'));
XStimParams.numitds = str2num(get(H.numitds,'String'));
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.curr_ILD = str2num(get(H.ILD,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.reset_flag = 0;

%Setup for ITDs used
itds  = round(XStimParams.loitd:...
   (XStimParams.hiitd - XStimParams.loitd)/(XStimParams.numitds-1):...
   XStimParams.hiitd);
%itds = XStimParams.numitds;

numreps = XStimParams.numreps;
numtrials = length(itds);

clear BUF

%Specify DAMA buffers
BUF.L1		= 1;
BUF.R1		= 2;
BUF.playseq_left	= 4;
BUF.playseq_right	= 5;
BUF.playspec	= 6;

%Make play sequence buffers
S232('allot16',BUF.playseq_left,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.L1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_left);

S232('allot16',BUF.playseq_right,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.R1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_right);


%Make play specification buffer
S232('allot16',BUF.playspec,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_left);
S232('make',1,BUF.playseq_right);
S232('make',2,0);
S232('qpop16',BUF.playspec);

%Check to see if earphone filters are being used
ephoneflag = 0;
if(exist1('H.fig') & XStimParams.ephone_flag == 1)
    ephoneflag = 1;
    ephonefname = [FN.ephone_path FN.ephone];
    ephonefilt_left  = (mtlrch(ephonefname,1))';
    ephonefilt_right = (mtlrch(ephonefname,2))';
end

%Check to see if location filters are being used
locflag = 0; str1 = [];
if(exist1('H.fig') & XStimParams.loc_flag == 1)
    locflag = 1;
    locfname = [FN.loc_path FN.loc];
    dirmat = sph2dbl(mtlrdir(locfname));
    locnum = find(dirmat(1,:) == XStimParams.loc_azel(2) & dirmat(2,:) == XStimParams.loc_azel(1));
    locfilt_left  = (mtlrch(locfname,2*locnum-1))';
    locfilt_right = (mtlrch(locfname,2*locnum))';
    str1 = ['hrtf: ' num2str(XStimParams.loc_azel(2)) 'el, ' num2str(XStimParams.loc_azel(1)) ' az'];
end

% write diary
if get(H.recorddata,'Value')
    tempstr = ['    ' str1 ' ITD2-test  atten: ' num2str(XStimParams.curr_ABL)];
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
if(ephoneflag == 0)
    S232('allot16',BUF.L1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
    S232('allot16',BUF.R1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
elseif(ephoneflag == 1 & locflag == 0)
    S232('allot16',BUF.L1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
    S232('allot16',BUF.R1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
elseif(ephoneflag == 1 & locflag == 1)
    S232('allot16',BUF.L1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.R1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
end

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
if(exist1('H.fig') & ephoneflag == 1 & locflag == 1)
    S232('PD1npts',1,(silence_len + DUR*(round(TDT.Fs/1000))) + TDT.itdfiltlen  + ...
        TDT.ephonefiltlen + TDT.hrtffiltlen);
elseif(exist1('H.fig') & ephoneflag == 1 & locflag == 0)
    S232('PD1npts',1,(silence_len + DUR*(round(TDT.Fs/1000))) + TDT.itdfiltlen  + TDT.ephonefiltlen);
elseif(exist1('H.fig') & ephoneflag == 0)
    S232('PD1npts',1,(silence_len + DUR*(round(TDT.Fs/1000))) + TDT.itdfiltlen);
end

%Load Earphone filters
if(exist1('H.itdfig') & get(H.ephoneuseit,'Value') &...
        get(H.locuseit,'Value')) %Use earphone and location filters
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
    S232('pushf',locfilt_left,length(locfilt_left));
    S232('PreLoadRaw',1,S232('DSPid',dspid_left_loc),'MONO','STACK','','',1.0,1.0,1);
    S232('dropall');
    %Load right
    S232('pushf',ephonefilt_right,length(ephonefilt_right));
    S232('PreLoadRaw',1,S232('DSPid',dspid_right_ephone),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    S232('dropall');
    S232('pushf',locfilt_right,length(locfilt_right));
    S232('PreLoadRaw',1,S232('DSPid',dspid_right_loc),'MONO','STACK','','',1.0,1.0,1);
    S232('dropall');
elseif(exist1('H.itdfig') & get(H.ephoneuseit,'Value') & get(H.locuseit,'Value')==0) %Use earphone filters only
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
    while (exist1('H.itdfig') & S232('APactive'))  pause(0); end
    S232('pushf',ephonefilt_left,length(ephonefilt_left));
    while (exist1('H.itdfig') & S232('APactive'))  pause(0); end
    S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    %Load right
    S232('pushf',ephonefilt_right,length(ephonefilt_right));
    while (exist1('H.itdfig') & S232('APactive'))  pause(0); end
    S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    while (exist1('H.itdfig') & S232('APactive'))  pause(0); end
elseif(exist1('H.itdfig') & get(H.ephoneuseit,'Value') == 0)
    S232('PD1clrsched',1);
    S232('PD1nstrms',1,2,0);
    S232('PD1specIB',1,S232('IB',0),S232('DAC',0));
    S232('PD1specIB',1,S232('IB',1),S232('DAC',1));
end

S232('PA4atten',1,abs(XStimParams.curr_ABL)-20);
S232('PA4atten',2,abs(XStimParams.curr_ABL)-20);


set(H.ephonefile,'Enable','off');
set(H.ephoneuseit,'Enable','off');
set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

%Set MII parameters
mii_us_per_sample = 10; %microsecond per sample
mii_separation = 100; %only take events separated by 100 samples (i.e., 1 ms)


ISI = XStimParams.test_ISI;
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000)); %correct for ITD filtlength
if(exist1('H.fig') & get(H.ephoneuseit,'Value')) %correct for ephone filtlength
    ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));
end

%Loop to make the stimuli we'll use
remreps = 1;
set(H.buildplay,'String','Building Stimuli');
set(H.remreps,'String',num2str(remreps));
repnum = 1;
buffcycle = 1;
finalspikematrix = zeros(1,length(itds));
numtrials = length(itds);
while (exist1('H.itdfig') & (repnum <= 1))
    %Randomize the stimuli
    randseq = randperm(numtrials);
    
    %Loop to make stimuli
    trialnum = 1;
    while (exist1('H.itdfig') & (trialnum <= numtrials))
        set(H.buildplay,'BackgroundColor','red');
        %Check for pause by user
        if pause_check  return; end
        
        switch get(H.stim_type,'Value')     
            case 1       % tone
                trial_left = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur);
                trial_right = trial_left;
            case 8         %Broadband Noise
                trial_left = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
                trial_right = trial_left;
            case 9         %Stimulus from file
                if(~exist1('stim_from_file'))
                    while exist1([FN.stim_path FN.stim]) ~=2
                        [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.*'],'Select stim#1 File');
                        if(FN.stim_path ~= 0)
                            set(H.stim_filename,'String',FN.stim);
                        end
                    end
                    fid = fopen([FN.stim_path FN.stim],'r');
                    stim_from_file = fread(fid,inf,'float');
                    fclose(fid);
                end
                trial_left = stim_from_file;
                trial_left = trial_left(:)';
                if(~exist1('stim_from_file2'))
                    while exist1([FN.stim_path2 FN.stim2]) ~=2
                        [FN.stim2,FN.stim_path2] = uigetfile([FN.stim_path2 '*.*'],'Select stim#2 File');
                        if(FN.stim_path2 ~= 0)
                            set(H.stim_filename2,'String',FN.stim2);
                        end
                    end
                    fid = fopen([FN.stim_path2 FN.stim2],'r');
                    stim_from_file2 = fread(fid,inf,'float');
                    fclose(fid);
                end
                trial_right = stim_from_file2;
                trial_right = trial_right(:)';
            otherwise
                set(H.stim_type,'Value',8);
                trial_left = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
                trial_right = trial_left;
                disp('Stimtype not supported for ITDtest2.  Reset for BROADBAND');    
        end
        
            % remove any DCoffset
    trial_left = trial_left - mom(trial_left,1);
    trial_right = trial_right - mom(trial_right,1);
    
    % normalize to ACPower (added 3/5/07 as in 2-source)
    if ~XStimParams.HiDynamicRange
        trial_left = trial_left / mom(trial_left,2);
        trial_right = trial_right / mom(trial_right,2);
    end

        % impose ILD NOT GETTING ILD by ATTENUATORS
        ILD = XStimParams.curr_ILD;
        % apply ILD half to each side
        trial_left  = trial_left * (10^(-ILD/(2*20)));
        trial_right = trial_right * (10^(ILD/(2*20)));
        
        % 3/7/07 multiply by freefield_Xfactors
        trial_left = trial_left * TDT.freefield_Lfactor;
        trial_right = trial_right * TDT.freefield_Rfactor;
        
        % remove any DCoffset
        trial_left = trial_left - mom(trial_left,1);
        trial_right = trial_right - mom(trial_right,1);
        
        % modulate stim1
        if get(H.stim_type,'Value') ~= 9 & ~strcmp(XStimParams.mod_type,'None')
            Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
            trial_right = trial_right .* Envelope(:)';
            trial_left = trial_left .* Envelope(:)';
        end
        
        %Ramp the stimuli
        ramp_time = 5; %ms
        [trial_left] = ramp_sound(trial_left,TDT.Fs,ramp_time);
        [trial_right] = ramp_sound(trial_right,TDT.Fs,ramp_time);
        
        %Add in the leading silent period
        trial_left =  [zeros(1,silence_len) trial_left];
        trial_right = [zeros(1,silence_len) trial_right];
        
        %Add in the trailing silent period
        trial_left =  [trial_left zeros(1,silence_len2)];
        trial_right = [trial_right zeros(1,silence_len2)];
        
        %Apply ITD filtering
        itdleft = 0; itdright = 0;
        ITD = round(itds(randseq(trialnum)));
        if(ITD < 0)
            itdleft = 0;
            itdright = abs(itds(randseq(trialnum)));
        elseif(ITD > 0)
            itdleft = abs(itds(randseq(trialnum)));
            itdright = 0;
        end
        eval(['load ' FN.itd_path 'itdfilt' num2str(itdleft)]);
        eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
        itd_filt_left = itd_filt_left/max(abs(itd_filt_left));
        eval(['load ' FN.itd_path 'itdfilt' num2str(itdright)]);
        eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
        itd_filt_right = itd_filt_right/max(abs(itd_filt_left));
        itdtrial_left = conv(trial_left,itd_filt_left);
        itdtrial_right = conv(trial_right,itd_filt_right);
        
        %pad for filters
        if(locflag == 1 & ephoneflag == 1)
            itdtrial_left = [itdtrial_left zeros(1,TDT.ephonefiltlen) zeros(1,TDT.hrtffiltlen)];
            itdtrial_right = [itdtrial_right zeros(1,TDT.ephonefiltlen) zeros(1,TDT.hrtffiltlen)];
        elseif(ephoneflag == 1)
            itdtrial_left = [itdtrial_left zeros(1,TDT.ephonefiltlen)];
            itdtrial_right = [itdtrial_right zeros(1,TDT.ephonefiltlen)];
        end
        
        % scale
        if 0
            ABAval = 0.5*(mom(itdtrial_left,2) + mom(itdtrial_right,2));
            scalefact = TDT.scalevalue/ABAval;
            itdtrial_left = round(scalefact*itdtrial_left * TDT.hrtf_scale);
            itdtrial_right = round(scalefact*itdtrial_right * TDT.hrtf_scale);
        end
        % 3/7/07
        itdtrial_left = itdtrial_left * TDT.scaleFactor;
        itdtrial_right = itdtrial_right * TDT.scaleFactor;
        
        
        %Fill stimulus buffers
        if(exist1('H.itdfig') & buffcycle == 1);
            S232('push16',itdtrial_left,length(itdtrial_left));
            S232('qpop16',BUF.L1);
            fname = ['L' num2str(itds(randseq(trialnum)))];
            evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
                    [' ''' stimuli_dir fname ''' '] ...
                    ',0);'];
            eval(evalstr);
            S232('push16',itdtrial_right,length(itdtrial_right));
            S232('qpop16',BUF.R1);
            fname = ['R' num2str(itds(randseq(trialnum)))];
            evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
                    [' ''' stimuli_dir fname ''' '] ...
                    ',0);'];
            eval(evalstr);
        end
        remtrials = numtrials - trialnum;
        set(H.remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        set(H.buildplay,'BackgroundColor','yellow');
        pause(0);
    end %end loop over trials
    
    remreps = 1;
    set(H.remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end %end loop over reps

%Put up raster plot
position = [1000 0 200 950];
xvals = 1;
yvals = itds;
xtext = '';
ytext = 'ITD (us)';
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

%Begin loop over experimental blocks (i.e., reps) and trials   
set(H.buildplay,'String','Playing Stimuli');
set(H.buildplay,'BackgroundColor','yellow');
set(H.remreps,'String',num2str(numreps));
repnum = 1;
datamatrix = [];

if(exist1('H.itdfig') & get(H.recorddata,'Value'))
    update_dataFN;
end
while (exist1('H.itdfig') & (repnum <= numreps))
    %Randomize the stimuli
    randseq = randperm(numtrials);   
    trialnum = 1;
    spikes_trial = [];
    itds_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    tic
    while (exist1('H.itdfig') & (trialnum <= numtrials))
        %Check for pause by user
        if pause_check  return; end
        
        fname = ['L' num2str(itds(randseq(trialnum)))];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' stimuli_dir fname ''' '] ');'];
        eval(evalstr);
        evalstr = ['S232(''qpop16'''  ',BUF.L1);'];
        eval(evalstr);
        fname = ['R' num2str(itds(randseq(trialnum)))];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' stimuli_dir fname ''' '] ');'];
        eval(evalstr);
        evalstr = ['S232(''qpop16'''  ',BUF.R1);'];
        eval(evalstr);
        
        S232('seqplay',BUF.playspec);
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
            
        while S232('PD1status',1) usec_delay(1000); end
        S232('PD1stop',1);
        
        m110dx( C_.STOP);
        spikes = m110dx( C_.DATA, 1000); %Take 100 spikes max
        ind = find(spikes ~= 0); %Get clock events that are spikes
        spikes = spikes(ind);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
        if(exist1('H.itdfig') & ~isempty(spikes)) 
            spikes = spikes(ind);
            spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
            itds_trial = [itds_trial;randseq(trialnum) * ones(size(spikes))];
            repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
            Nspikes = [Nspikes; ones(length(spikes),1)*length(spikes)];
        end
        
        
        finalspikematrix(randseq(trialnum)) = ...
            finalspikematrix(randseq(trialnum)) + ...
            length(spikes);
        
        if pause_check    return; end
        
        remtrials = numtrials - trialnum;
        set(H.remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end %end loop over trials
    
    %Plot all the spikes from trials
    plotraster(H.rasterfig,...
        spikes_trial,...
        ones(size(itds_trial)),...
        itds_trial,...
        repnum,...
        tot_dur,...
        stim_dur,...
        numreps,...
        xoffset,...
        yoffset,...
        row_space,...
        col_space);
    
    %Record Data
    if(exist1('H.itdfig') & get(H.recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial itds_trial]];
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
plot(itds,finalspikematrix,'g-*','LineWidth',1.5);
set(H.axes,'Color','black');
xlabel('ITD (us)'); ylabel('Spike Rate (spikes/stim)');
title('ITD Test');

set(H.buildplay,'String','Build/Play status');
set(H.exitITDtest2,'Visible','on');
set(H.resetITDtest2,'Visible','on');

% increment testnum
if(exist1('H.itdfig') & get(H.recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
end
%%%%%%%%%
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

global XStimParams
global H

flag = 0;
%Check for pause by user
while (exist1('H.itdfig') & get(H.pauseITDtest2,'Value'))
    pause(0);
    if(~exist1('H.itdfig')) return; end         
    set(H.exitITDtest2,'Visible','on');
    set(H.resetITDtest2,'Visible','on');
    if(exist1('H.itdfig') & get(H.resetITDtest2,'Value') == 1)
        set(H.resetITDtest2,'Value',0);
        set(H.pauseITDtest2,'Value',0);
        Reset_ITDtest2;  flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end