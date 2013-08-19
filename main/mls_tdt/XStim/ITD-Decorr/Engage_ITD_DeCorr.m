function [] = Engage_itd_decorr()

global H
global XStimParams
global TDT
global FN
global C_
global M_

%Engage_itd_decorr

%*******************************************************************************
%	The itd_decorr operation
% this version allows different sounds for each ear
% will play interleaved trials of different filenames
% chosen from the matrix of stimuli
% NOTE: testtypes in DATA file assigned by XStimParams.stims_to_play down column 1, down column2,...
% so that e.g. if three highest priority FNs played, 
% 1)Lo (plotted blue)  2) Hi (green)  3) Hi & Lo (cyan)
%*******************************************************************************
stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);
disp('This is an ITD_decorr test')
marker_color= 'bgcmywbgcmyw';

% be sure stimFNs picked
dur = [XStimParams.curr_stimdur XStimParams.curr_stimdur];
if any(XStimParams.stims_to_play(:,2)) & (isempty(FN.stim4L) | isempty(FN.stim4R))
    set(H.stim_FN4pb,'value',1);
    setinfo_itd_DeCorr;
end
if any(XStimParams.stims_to_play(:,3)) & (isempty(FN.stim5L) | isempty(FN.stim5R))
    set(H.stim_FN5pb,'value',1);
    setinfo_itd_DeCorr;
end
if any(XStimParams.stims_to_play(:,4)) & (isempty(FN.stim6L) | isempty(FN.stim6R))
    set(H.stim_FN6pb,'value',1);
    setinfo_itd_DeCorr;
end
if any(XStimParams.stims_to_play(2,:)) & (isempty(FN.stim7L) | isempty(FN.stim7R))
    set(H.stim_FN7pb,'value',1);
    setinfo_itd_DeCorr;
end
if any(XStimParams.stims_to_play(3,:)) & (isempty(FN.stim8L) | isempty(FN.stim8R))
    set(H.stim_FN8pb,'value',1);
    setinfo_itd_DeCorr;
end
if any(XStimParams.stims_to_play(4,:)) & (isempty(FN.stim9L) | isempty(FN.stim9R))
    set(H.stim_FN9pb,'value',1);
    setinfo_itd_DeCorr;
end

%Put parameters into XStimParams
XStimParams.test_type = 'ITD_decorr';
XStimParams.stim_type = 'File';
XStimParams.loitd = str2num(get(H.lowitd,'String'));
XStimParams.hiitd = str2num(get(H.highitd,'String'));
XStimParams.numitds = str2num(get(H.numitds,'String'));
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.curr_ILD = str2num(get(H.ILD,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.reset_flag = 0;

%Setup for ABLs used
itds  = round(XStimParams.loitd: ...
   (XStimParams.hiitd - XStimParams.loitd)/(XStimParams.numitds-1): ...
   XStimParams.hiitd);
nITDs = length(itds);
%itds = XStimParams.numitds;

clear BUF
%Specify DAMA buffers
BUF.stimleft_1		= 1;
BUF.stimright_1		= 2;
BUF.playseq_left_1	= 4;
BUF.playseq_right_1	= 5;
BUF.playspecbuf_1	= 6;
BUF.stimleft_2		= 7;
BUF.stimright_2		= 8;
BUF.playseq_left_2	= 9;
BUF.playseq_right_2	= 10;
BUF.playspecbuf_2	= 11;

%Make play sequence buffers
S232('allot16',BUF.playseq_left_1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.stimleft_1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_left_1);

S232('allot16',BUF.playseq_right_1,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.stimright_1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_right_1);

S232('allot16',BUF.playseq_left_2,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.stimleft_2);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_left_2);

S232('allot16',BUF.playseq_right_2,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.stimright_2);
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
    str1 = ['hrtf: ' num2str(XStimParams.loc_azel(2)) 'el, ' num2str(XStimParams.loc_azel(1)) ' az'];
end

% write diary
if get(H.recorddata,'Value')
    tempstr = ['    ' str1 ' ITD-decorr-test  atten: ' num2str(XStimParams.curr_ABL)];
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
    S232('allot16',BUF.stimleft_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
    S232('allot16',BUF.stimright_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
    S232('allot16',BUF.stimleft_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
    S232('allot16',BUF.stimright_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
elseif(ephoneflag == 1 & locflag == 0)
    S232('allot16',BUF.stimleft_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
    S232('allot16',BUF.stimright_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
    S232('allot16',BUF.stimleft_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
    S232('allot16',BUF.stimright_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
elseif(ephoneflag == 1 & locflag == 1)
    S232('allot16',BUF.stimleft_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.stimright_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.stimleft_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.stimright_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
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
if(exist1('H.itd_decorrfig') & get(H.ephoneuseit,'Value') &...
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
elseif(exist1('H.itd_decorrfig') & get(H.ephoneuseit,'Value') & get(H.locuseit,'Value')==0) %Use earphone filters only
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
    while (exist1('H.itd_decorrfig') & S232('APactive'))  pause(0); end
    S232('pushf',ephonefilt_left,length(ephonefilt_left));
    while (exist1('H.itd_decorrfig') & S232('APactive'))  pause(0); end
    S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    %Load right
    S232('pushf',ephonefilt_right,length(ephonefilt_right));
    while (exist1('H.itd_decorrfig') & S232('APactive'))  pause(0); end
    S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    while (exist1('H.itd_decorrfig') & S232('APactive'))  pause(0); end
elseif(exist1('H.itd_decorrfig') & get(H.ephoneuseit,'Value') == 0)
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
nTests = length(find(XStimParams.stims_to_play));
remtrials = nITDs * nTests;
set(H.buildplay,'String','Building Stimuli');
set(H.remreps,'String',num2str(1));
set(H.remtrials,'string',num2str(remtrials));

trialnum = 0;
set(H.buildplay,'BackgroundColor','red');
for iHi = 1:4             % loop through the HiFrequency FNs
    for iLo = 1:4         % Loop through the LoFrequency FNs
        if XStimParams.stims_to_play(iLo,iHi)
            
            %Check for pause by user
            if pause_check  return; end
            
            % Left ear
            % read HiFreq stimfile first
            if iHi > 1
                temp_FN = [FN.stim_path eval(['FN.stim'  num2str(iHi+2) 'L'])];
                fid = fopen(temp_FN,'r');
                trial_left = fread(fid,inf,'float');
                fclose(fid);
                trial_left = trial_left(:)';
                trial_left = trial_left - mom(trial_left,1);
            else
                trial_left = zeros(DUR*(round(TDT.Fs/1000)),1);
            end
            % read LoFreq stimfile
            if iLo > 1
                temp_FN = [FN.stim_path eval(['FN.stim'  num2str(iLo+5) 'L'])];
                fid = fopen(temp_FN,'r');
                stim_from_file = fread(fid,inf,'float');
                fclose(fid);
                trial_left = trial_left + stim_from_file;
            end
            % Right ear
            % read HiFreq stimfile first
            if iHi > 1
                temp_FN = [FN.stim_path eval(['FN.stim'  num2str(iHi+2) 'R'])];
                fid = fopen(temp_FN,'r');
                trial_right = fread(fid,inf,'float');
                fclose(fid);
                trial_right = trial_right(:)';
                trial_right = trial_right - mom(trial_right,1);
            else
                trial_right = zeros(DUR*(round(TDT.Fs/1000)),1);     %%%%%?
            end
            % read LoFreq stimfile
            if iLo > 1
                temp_FN = [FN.stim_path eval(['FN.stim'  num2str(iLo+5) 'R'])];
                fid = fopen(temp_FN,'r');
                stim_from_file = fread(fid,inf,'float');
                fclose(fid);
                trial_right = trial_right + stim_from_file;
            end
 
            % scale to ACpower == 1;
            trial_left = trial_left / mom(trial_left,2);
            trial_right = trial_right / mom(trial_right,2);
            
            ILD = XStimParams.curr_ILD;
            % apply ILD half to each side
            trial_left  = stim * (10^(-ILD/(2*20)));
            trial_right = stim * (10^(ILD/(2*20)));
            
            % remove any DCoffset
            trial_left = trial_left - mom(trial_left,1);
            trial_right = trial_right - mom(trial_right,1);
            
            % 3/7/07 multiply by freefield_Xfactors
            trial_left = trial_left * TDT.freefield_Lfactor;
            trial_right = trial_right * TDT.freefield_Rfactor;
            
            %Ramp the stimuli
            ramp_time = 5; %ms
            [trial_left] = ramp_sound(trial_left,TDT.Fs,ramp_time);
            [trial_right] = ramp_sound(trial_right,TDT.Fs,ramp_time);
            
            %Add in the leading silent period
            trial_left =  [zeros(1,silence_len) trial_left(:)'];
            trial_right = [zeros(1,silence_len) trial_right(:)'];
            
            %Add in the trailing silent period
            trial_left =  [trial_left zeros(1,silence_len2)];
            trial_right = [trial_right zeros(1,silence_len2)];
            
            %pad for filters
            if(locflag == 1 & ephoneflag == 1)
                trial_left = [trial_left zeros(1,TDT.ephonefiltlen) zeros(1,TDT.hrtffiltlen)];
                trial_right = [trial_right zeros(1,TDT.ephonefiltlen) zeros(1,TDT.hrtffiltlen)];
            elseif(ephoneflag == 1)
                trial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
                trial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
            end
            
            %%%%%%%%%Apply ITD filtering
            for iITD = 1:nITDs
                trialnum = trialnum + 1;        % over all stim combinations
                itdleft = 0; itdright = 0;
                ITD = round(itds(iITD));
                if(ITD < 0)
                    itdright = abs(itds(iITD));
                elseif(ITD > 0)
                    itdleft = abs(itds(iITD));
                end
                eval(['load ' FN.itd_path 'itdfilt' num2str(itdleft)]);
                eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
                itd_filt_left = itd_filt_left/max(abs(itd_filt_left));
                eval(['load ' FN.itd_path 'itdfilt' num2str(itdright)]);
                eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
                itd_filt_right = itd_filt_right/max(abs(itd_filt_left));
                itdtrial_left = conv(trial_left,itd_filt_left);
                itdtrial_right = conv(trial_right,itd_filt_right);
                
                % remove DCoffset
                itdtrial_left = itdtrial_left - round(mean(itdtrial_left));
                itdtrial_right = itdtrial_right - round(mean(itdtrial_right));
                
                % 3/7/07
                itdtrial_left = itdtrial_left * TDT.scaleFactor;
                itdtrial_right = itdtrial_right * TDT.scaleFactor;
                
                %Fill stimulus buffers
                S232('push16',itdtrial_left,length(itdtrial_left));
                S232('qpop16',BUF.stimleft_1);
                fname = ['stimleft_' num2str(iHi) '_' num2str(iLo) '_' num2str(itds(iITD))];
                evalstr = ['S232(''dama2disk16'',BUF.stimleft_1,' ...
                        [' ''' stimuli_dir fname ''' '] ',0);'];
                eval(evalstr);
                FN_listL{trialnum} = fname;
                
                S232('push16',itdtrial_right,length(itdtrial_right));
                S232('qpop16',BUF.stimright_1);
                fname = ['stimright_' num2str(iHi) '_' num2str(iLo) '_' num2str(itds(iITD))];
                evalstr = ['S232(''dama2disk16'',BUF.stimright_1,' ...
                        [' ''' stimuli_dir fname ''' '] ',0);'];
                eval(evalstr);
                FN_listR{trialnum} = fname;
                remtrials = remtrials-1;
                set(H.remtrials,'string',num2str(remtrials));
                pause(0);
            end % iITD = 1:length(itds)
            pause(0);
        end % XStimParams.stims_to_play(iLo,iHi)
    end         % iLo
end             % iHi
numtrials = trialnum;
%%%% finished making stimuli
numreps = XStimParams.numreps;
set(H.buildplay,'BackgroundColor','yellow');

%Put up raster plot
position = [1000 50 200 900];
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
    xvals, yvals,...
    xtext, ytext,...
    xoffset, yoffset,...
    numreps,...
    tot_dur, stim_dur,...
    startcount, endcount,...
    row_space, col_space);

% for plotting
finalspikematrix = zeros(nTests,length(itds));                   %%%%%??? more of these??

%Begin loop over experimental blocks (i.e., reps) and trials   
set(H.buildplay,'String','Playing Stimuli');
set(H.buildplay,'BackgroundColor','yellow');
set(H.remreps,'String',num2str(numreps));
repnum = 1;
buffcycle = 1;
datamatrix = [];

if(exist1('H.itd_decorrfig') & get(H.recorddata,'Value'))
    update_dataFN;
end
while (exist1('H.itd_decorrfig') & (repnum <= numreps))
    %Randomize the stimuli
    randseq = randperm(numtrials);   
    trialnum = 1;
    spikes_trial = [];
    itds_trial = [];
    repnum_trial = [];
    Nspikes = [];
    iHi_trial = [];
    iLo_trial = [];
    test_type_trial = [];
    while (exist1('H.itd_decorrfig') & (trialnum <= numtrials+1))
        %Check for pause by user
        if pause_check  return; end
        
        if(trialnum <= numtrials)
            iITD = randseq(trialnum);         % the randomly picked FN# to be loaded
            if(buffcycle == 1);
                fname = FN_listL{iITD};
                evalstr = ['S232(''disk2dama16'',BUF.stimleft_1,' ...
                        [' ''' stimuli_dir fname ''' '] ',0);'];
                eval(evalstr);
                fname = FN_listR{iITD};
                evalstr = ['S232(''disk2dama16'',BUF.stimright_1,' ...
                        [' ''' stimuli_dir fname ''' '] ',0);'];
                eval(evalstr);
            elseif(buffcycle == 2);
                fname = FN_listL{iITD};
                evalstr = ['S232(''disk2dama16'',BUF.stimleft_2,' ...
                        [' ''' stimuli_dir fname ''' '] ',0);'];
                eval(evalstr);
                fname = FN_listR{iITD};
                evalstr = ['S232(''disk2dama16'',BUF.stimright_2,' ...
                        [' ''' stimuli_dir fname ''' '] ',0);'];
                eval(evalstr);
            end
        end
        
        %Wait till PD1 is finished
        %while S232('PD1status',1) usec_delay(1000); end
        usec_delay((DUR+100)*1000);
        S232('PD1stop',1);
        %Check for pause by user
        if pause_check  return; end
        
        %Stop the m110 and get spikes
        if(trialnum > 1)
            iITD = randseq(trialnum-1);        % the randomly picked FN# that was just played
            iHi = str2num(FN_listL{iITD}(10));
            iLo = str2num(FN_listL{iITD}(12));
            itd = itds(mod(iITD,nITDs)+1);

            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, 1000); %Take 100 spikes max
            spikes = spikes(find(spikes ~= 0)); %Get clock events that are spikes
            ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
            testtype = find(find(XStimParams.stims_to_play)==sub2ind([4 4],iLo,iHi));
            if(exist1('H.itd_decorrfig') & ~isempty(spikes)) 
                spikes = spikes(ind);
                spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
                itds_trial = [itds_trial;(mod(iITD,nITDs)+1) * ones(size(spikes))];         %%%% ????
                repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
                Nspikes = [Nspikes; ones(length(spikes),1)*length(spikes)];
                iHi_trial = [iHi_trial; iHi * ones(size(spikes))];
                iLo_trial = [iLo_trial; iLo * ones(size(spikes))];
                test_type_trial = [test_type_trial; testtype * ones(size(spikes))];
            end
        end
        
        if(exist1('H.itd_decorrfig') & trialnum <= numtrials)
            
            if(exist1('H.itd_decorrfig') & buffcycle == 1)
                S232('seqplay',BUF.playspecbuf_1);
            elseif(exist1('H.itd_decorrfig') & buffcycle == 2)
                S232('seqplay',BUF.playspecbuf_2);
            end
            S232('PD1arm',1);
            
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
            testtype = find(find(XStimParams.stims_to_play)==sub2ind([4 4],iLo,iHi));
            ind0 = find(spikes < 100*DUR_silence & spikes > 0);            % spikes preceding sound onset
            ind1 = find(spikes > 100*DUR_silence & spikes <= (DUR_silence + DUR)*100);   % spikes during sound
            evoked =  length(ind1) - length(ind0);
            finalspikematrix(testtype,mod(iITD,nITDs)+1) = ...
                finalspikematrix(testtype,mod(iITD,nITDs)+1) + evoked;
        end
        
        remtrials = numtrials - trialnum;
        set(H.remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end %end loop over trials
    
    %Plot all the spikes from trials
    for itest = 1:nTests
        ind = find(test_type_trial == itest);
    plotraster(H.rasterfig,...
        spikes_trial(ind),...
        ones(size(itds_trial(ind))),...
        itds_trial(ind),...
        repnum,...
        tot_dur, stim_dur,...
        numreps,...
        xoffset, yoffset,...
        row_space, col_space, marker_color(itest));
    end
    
    %Record Data
    if(exist1('H.itd_decorrfig') & get(H.recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial itds_trial iHi_trial iLo_trial]];
        record_data3(XStimParams,datamatrix);
    end
    
    
    %Plot final spike rate figure
    if exist1('H.finalspikeratefig')
        close(H.finalspikeratefig);
    end
    H.finalspikeratefig = figure('Position',[480 300 500 450],...
        'Name','Current Test Results',...
        'NumberTitle','off');
    H.axes = axes;
    hold on
    for iplot = 1:nTests
        plot(itds,finalspikematrix(iplot,:)/repnum,'-*','color',marker_color(iplot),'LineWidth',1.5);
    end
    set(H.axes,'Color','black');
    xlabel('ITD (us)'); ylabel('Spike Rate (spikes/stim)');
    title('ITD decorr Test');
    axis([min(itds) max(itds) min1([0 floor(min1(finalspikematrix/repnum))]) max1([1 ceil(max1(finalspikematrix/repnum))]) ])
    legh = legend(XStimParams.legendstr); 
    set(legh,'position',[.6 (1 - length(XStimParams.legendstr)*.04) .4 length(XStimParams.legendstr)*.04]);
    set(legh,'color',[.9 .9 .9]);
    %%%%%

    set(H.remreps,'String',num2str(numreps - repnum));
    repnum = repnum + 1;
    pause(0);
end %end loop over reps


set(H.buildplay,'String','Build/Play status');
set(H.exititd_decorr,'Visible','on');
set(H.resetitd_decorr,'Visible','on');

% increment testnum
if(exist1('H.itd_decorrfig') & get(H.recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flag] = pause_check

global XStimParams
global H

flag = 0;
%Check for pause by user
while (exist1('H.itd_decorrfig') & get(H.pauseitd_decorr,'Value'))
    pause(0);
    if(~exist1('H.itd_decorrfig')) return; end         
    set(H.exititd_decorr,'Visible','on');
    set(H.resetitd_decorr,'Visible','on');
    if(exist1('H.itd_decorrfig') & get(H.resetitd_decorr,'Value') == 1)
        set(H.resetitd_decorr,'Value',0);
        set(H.pauseitd_decorr,'Value',0);
        Reset_ITD_decorr;  flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end