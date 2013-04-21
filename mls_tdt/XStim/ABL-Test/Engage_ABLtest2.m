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

%Put parameters into XStimParams
test_val = get(H.test_type,'Value');
test_type = get(H.test_type,'String');
XStimParams.test_type = deblank(test_type(test_val,:));
clear test_type
XStimParams.loabl = str2num(get(H.lowabl,'String'));
XStimParams.hiabl = str2num(get(H.highabl,'String'));
XStimParams.numabls = str2num(get(H.numabls,'String'));
XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
XStimParams.curr_ILD = str2num(get(H.ILD,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.reset_flag = 0;

%Setup for ABLs used
abls  = round(XStimParams.loabl:...
    (XStimParams.hiabl - XStimParams.loabl)/(XStimParams.numabls-1):...
    XStimParams.hiabl);

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

%Check to see if earphone filters are being used
ephoneflag = 0;
if(exist1('H.fig') & XStimParams.ephone_flag == 1)
    ephoneflag = 1;
    ephonefname = [FN.ephone_path FN.ephone];
    ephonefilt_left  = (mtlrch(ephonefname,1))';
    ephonefilt_right = (mtlrch(ephonefname,2))';
end

%Check to see if location filters are being used
locflag = 0;
if(exist1('H.fig') & XStimParams.loc_flag == 1)
    locflag = 1;
    locfname = [FN.loc_path FN.loc];
    dirmat = sph2dbl(mtlrdir(locfname));
    locnum = find(dirmat(1,:) == XStimParams.loc_azel(2) & dirmat(2,:) == XStimParams.loc_azel(1));
    locfilt_left  = (mtlrch(locfname,2*locnum-1))';
    locfilt_right = (mtlrch(locfname,2*locnum))';
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
    S232('allotf',BUF.left_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
    S232('allotf',BUF.right_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
    S232('allotf',BUF.left_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
    S232('allotf',BUF.right_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
elseif(ephoneflag == 1 & locflag == 0)
    S232('allotf',BUF.left_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
    S232('allotf',BUF.right_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
    S232('allotf',BUF.left_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
    S232('allotf',BUF.right_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
elseif(ephoneflag == 1 & locflag == 1)
    S232('allotf',BUF.left_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allotf',BUF.right_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allotf',BUF.left_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allotf',BUF.right_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
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
if(exist1('H.ablfig') & get(H.ephoneuseit,'Value') &...
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
elseif(exist1('H.ablfig') & get(H.ephoneuseit,'Value') & get(H.locuseit,'Value')==0) %Use earphone filters only
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
    while (exist1('H.ablfig') & S232('APactive'))  pause(0); end
    S232('pushf',ephonefilt_left,length(ephonefilt_left));
    while (exist1('H.ablfig') & S232('APactive'))  pause(0); end
    S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    %Load right
    S232('pushf',ephonefilt_right,length(ephonefilt_right));
    while (exist1('H.ablfig') & S232('APactive'))  pause(0); end
    S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    while (exist1('H.ablfig') & S232('APactive'))  pause(0); end
elseif(exist1('H.ablfig') & get(H.ephoneuseit,'Value') == 0)
    S232('PD1clrsched',1);
    S232('PD1nstrms',1,2,0);
    S232('PD1specIB',1,S232('IB',0),S232('DAC',0));
    S232('PD1specIB',1,S232('IB',1),S232('DAC',1));
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
if(exist1('H.fig') & get(H.ephoneuseit,'Value')) %correct for ephone filtlength
    ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));
end

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
    %Lright - Lleft = 20*log10(Ampright/Ampleft), Hartmann p. 29
    ILD = XStimParams.curr_ILD;
    mult_fact = 10^(ILD/20);
    [bbnoise] = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
    trial_left = bbnoise/sqrt(mult_fact);
    trial_right = sqrt(mult_fact) * bbnoise;
    
    % remove any DCoffset
    trial_left = trial_left - mom(trial_left,1);
    trial_right = trial_right - mom(trial_right,1);
    
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
    abltrial_left = conv(trial_left,itd_filt_left);
    abltrial_right = conv(trial_right,itd_filt_right);
    
    %Correct for filters
    if(locflag == 1 & ephoneflag == 1)
        abltrial_left = [abltrial_left zeros(1,TDT.ephonefiltlen) zeros(1,hrtffiltlen)];
        abltrial_right = [abltrial_right zeros(1,TDT.ephonefiltlen) zeros(1,hrtffiltlen)];
    elseif(ephoneflag == 1)
        abltrial_left = [abltrial_left zeros(1,TDT.ephonefiltlen)];
        abltrial_right = [abltrial_right zeros(1,TDT.ephonefiltlen)];
    end
    
    %%%%%   scaling moved here on 3/17/03
    %ABAval = 0.5*(max(abs(abltrial_left)) + max(abs(abltrial_right)));
    abltrial_left = abltrial_left - round(mean(abltrial_left));
    abltrial_right = abltrial_right - round(mean(abltrial_right));
    ABAval = 0.5*(mom(abltrial_left,2) + mom(abltrial_right,2));
    scalefact = TDT.scalevalue/ABAval;
    abltrial_left = round(scalefact*abltrial_left * TDT.hrtf_scale);
    abltrial_right = round(scalefact*abltrial_right * TDT.hrtf_scale);
    
    
    % add in ABL        Aug 30, 2006    
    factor = 10^((abs(abls(trialnum)) + max(abls))/20);     % max(abls) is negative
    
    trial_left = abltrial_left/factor;
    trial_right = abltrial_right/factor;
    
    % write files
    S232('pushf',trial_left,length(trial_left));
    %S232('qpopf',BUF.left_1);
    fname = ['ABL.left_' num2str(abs(abls(trialnum)))];
    evalstr = ['S232(''popdiskf'',' ...
            [' ''' FN.temp_stim_path fname ''' '] ');'];
    eval(evalstr);
    
    S232('pushf',trial_right,length(trial_right));
    %S232('qpopf',BUF.right_1);
    fname = ['ABL.right_' num2str(abs(abls(trialnum)))];
    
    evalstr = ['S232(''popdiskf'',' ...
            [' ''' FN.temp_stim_path fname ''' '] ');'];
    eval(evalstr);
    
    
    %evalstr = ['S232(''dama2disk16'',BUF.right_1,' ...
    %        [' ''' FN.temp_stim_path fname ''' '] ',0);'];
    %eval(evalstr);
    
    %Plot PSD
    if(exist1('H.ablfig') & get(H.plotpsd,'Value') == 1)
    elseif(exist1('H.ablfig') & exist1('H.pl') & get(H.plotpsd,'Value') == 0)
        set(H.psdaxes,'Visible','off');
        set(H.pl,'Visible','off');
        set(H.pr,'Visible','off');
    end
    set(H.remtrials,'String',num2str(numtrials-trialnum));
    trialnum = trialnum + 1;
    pause(0);
    
end %end loop over trials to build stimuli

S232('PA4atten',1,min(abs(abls))-20);
S232('PA4atten',2,min(abs(abls))-20);


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
                    fname = ['ABL.left_' num2str(abs(abls(randseq(trialnum))))];
                    evalstr = ['disktodamaf(' ...
                            [' ''' FN.temp_stim_path fname ''' ']  ',BUF.left_1,0,[],[]);'];
                    eval(evalstr);
                    fname = ['ABL.right_' num2str(abs(abls(randseq(trialnum))))];
                    evalstr = ['disktodamaf(' ...
                            [' ''' FN.temp_stim_path fname ''' ']  ',BUF.right_1,0,[],[]);'];
                    eval(evalstr);
                elseif(buffcycle == 2);
                    fname = ['ABL.left_' num2str(abs(abls(randseq(trialnum))))];
                    evalstr = ['disktodamaf(' ...
                            [' ''' FN.temp_stim_path fname ''' ']  ',BUF.left_2,0,[],[]);'];
                    eval(evalstr);
                    fname = ['ABL.right_' num2str(abs(abls(randseq(trialnum))))];
                    evalstr = ['disktodamaf(' ...
                            [' ''' FN.temp_stim_path fname ''' ']  ',BUF.right_2,0,[],[]);'];
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

S232('PA4atten',1,80);
S232('PA4atten',2,80);

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