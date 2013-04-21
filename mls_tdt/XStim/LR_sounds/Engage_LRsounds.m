function [] = Engage_LRsounds()

global H
global XStimParams
global TDT
global FN
global C_
global M_

%Engage_LRsounds

%*******************************************************************************
% The LRsounds operation allows different sounds for each ear
% taken from files
%*******************************************************************************
stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);
disp('This is an LRsounds test')

% be sure stimFNs picked
if isempty(FN.stim)
    set(H.LR_stimFNpb(1),'value',1);
    setinfo_LRsounds;
end
if isempty(FN.stim2)
    set(H.LR_stimFNpb(2),'value',1);
    setinfo_LRsounds;
end
fid = fopen([FN.stim_path FN.stim],'r');
stim1 = fread(fid,inf,'float');
fclose(fid);
temp = mom(stim1,1);
if abs(temp) >.01
    B = questdlg(['The LEFT stimulus is offset from zero by ' num2str(temp), ', do you want to scale it to zero-offset?'], ...
        'Stimulus Offset', 'rescale','leave as is','rescale');
    if strcmp(B,'rescale')
        stim1 = stim1 - temp;
    end
end
temp = mom(stim1,2);
if temp > 1
    B = questdlg(['The LEFT stimulus has ACpower of ' num2str(temp), ', do you want to scale it to ACpower==1?'], ...
        'High ACpower', 'rescale','leave as is','rescale');
    if strcmp(B,'rescale')
        stim1 = stim1/temp;
    end
end

fid = fopen([FN.stim_path2 FN.stim2],'r');
stim2 = fread(fid,inf,'float');
fclose(fid);
temp = mom(stim2,1);
if abs(temp) >.01
    B = questdlg(['The RIGHT stimulus is offset from zero by ' num2str(temp), ', do you want to scale it to zero-offset?'], ...
        'Stimulus Offset', 'rescale','leave as is','rescale');
    if strcmp(B,'rescale')
        stim2 = stim2 - temp;
    end
end
temp = mom(stim2,2);
if temp > 1
    B = questdlg(['The RIGHT stimulus has ACpower of ' num2str(temp), ', do you want to scale it to ACpower==1?'], ...
        'High ACpower', 'rescale','leave as is','rescale');
    if strcmp(B,'rescale')
        stim2 = stim2/temp;
    end
end

%Put parameters into XStimParams
XStimParams.test_type = 'LRsounds';
XStimParams.stim_type = 'File';
XStimParams.reset_flag = 0;
nITDs = length(XStimParams.itds);
nILDs = length(XStimParams.ilds);
nABLs = length(XStimParams.abls);

% load HRTFs
ephoneflag = 0;
if exist1('GUI2.locations')
    XStimParams.locations = GUI2.locations;
    nLocs = size(GUI2.locations,2);
    HRTFfiletype = testHRTFfiletype(GUI.hrtfPATH, GUI.hrtfFN);
    if HRTFfiletype == 1
        HRTFfname = [GUI.hrtfPATH, GUI.hrtfFN];
        HRTFdir = mtlrdir(HRTFfname);
        for iLoc = 1:nLocs
            ind = find(HRTFdir(1,:)==XStimParams.locations(1,iLoc) & HRTFdir(2,:)==XStimParams.locations(2,iLoc));
            hrTF1(iLoc,:) = mtlrch(HRTFfname,ind*2-11);
            hrTF2(iLoc,:) = mtlrch(HRTFfname,ind*2);
        end
    else
        dir = [];
        load [GUI.hrtfPATH, GUI.hrtfFN]
        HRTFdir = dir;
        for iLoc = 1:nLocs
            ind = find(HRTFdir(1,:)==XStimParams.locations(1,iLoc) & HRTFdir(2,:)==XStimParams.locations(2,iLoc));
            hrTF1(iLoc,:) = TF1(ind,:);
            hrTF2(iLoc,:) = TF2(ind,:);
        end
    end
    clear dir TF*
    
    if findstr(GUI2.hrtfFN,'std') | XStimParams.ephone_flag
        ephoneflag = 1;
        EPHONEfiletype = testHRTFfiletype(FN.ephone_path, FN.ephone);
        if EPHONEfiletype ==1
            ephonefname = [FN.ephone_path FN.ephone];
            ephonefilt_left  = (mtlrch(ephonefname,1))';
            ephonefilt_right = (mtlrch(ephonefname,2))';
        else
            dir = [];
            load [FN.ephone_path FN.ephone]
            ephonefilt_left  = TF1;
            ephonefilt_right = TF2;
            clear TF* dir
        end
        clear EPHONEfiletype
    end
else
    XStimParams.locations = [];
    nLocs = 0;
end

clear BUF
%Specify DAMA buffers
BUF.L1		        = 1;
BUF.R1		        = 2;
BUF.playseq_left	= 4;
BUF.playseq_right	= 5;
BUF.playspecbuf	    = 6;

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
S232('allot16',BUF.playspecbuf,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_left);
S232('make',1,BUF.playseq_right);
S232('make',2,0);
S232('qpop16',BUF.playspecbuf);

% silence prior to stimulus
silence_lead_pts = (XStimParams.silence_lead * round(TDT.Fs/1000));
% silence after stimulus 
silence_trail_pts = (XStimParams.silence_trail * round(TDT.Fs/1000));

%Make Stimulus buffers
DUR = str2num(get(H.DUR,'String'));
if(ephoneflag == 0)
    Npts = silence_lead_pts + DUR*round(TDT.Fs/1000) + silence_trail_pts + TDT.itdfiltlen;
elseif(ephoneflag == 1 & locflag == 0)
    Npts = silence_lead_pts + DUR*round(TDT.Fs/1000) + silence_trail_pts + TDT.itdfiltlen + TDT.ephonefiltlen;
elseif(ephoneflag == 1 & locflag == 1)
    Npts = silence_lead_pts + DUR*round(TDT.Fs/1000) + silence_trail_pts + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen;
end
S232('allot16', BUF.L1, Npts);
S232('allot16', BUF.R1, Npts);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,Npts);

%Load Earphone filters
if (ephoneflag & nLocs)         % Use earphone and HRTFs
    S232('PD1clrsched',1);
    S232('PD1nstrms',1,2,0);
    S232('PD1resetDSP',1,hex2dec('FFF'));
    S232('dropall');
    %Make connections for left ear
    S232('PD1addsimp',1,S232('DSPout',2),S232('DAC',0));            %ephone to DAC0
    S232('PD1addsimp',1,S232('DSPout',0),S232('DSPin',2));          %loc to ephone
    S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left_loc));  %IB to loc
    %Make connections for right ear
    S232('PD1addsimp',1,S232('DSPout',3),S232('DAC',1));            %ephone to DAC1
    S232('PD1addsimp',1,S232('DSPout',1),S232('DSPin',dspid_right_ephone)); %loc to ephone
    S232('PD1specIB',1,S232('IB',1),S232('DSPin',1));               %IB to loc
    
    %Load left      
    S232('pushf',ephonefilt_left,length(ephonefilt_left));
    S232('PreLoadRaw',1,S232('DSPid',2),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    S232('dropall');
    %S232('pushf',locfilt_left,length(locfilt_left));
    %S232('PreLoadRaw',1,S232('DSPid',0),'MONO','STACK','','',1.0,1.0,1);
    %S232('dropall');
    %Load right
    S232('pushf',ephonefilt_right,length(ephonefilt_right));
    S232('PreLoadRaw',1,S232('DSPid',3),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    S232('dropall');
    %S232('pushf',locfilt_right,length(locfilt_right));
    %S232('PreLoadRaw',1,S232('DSPid',1),'MONO','STACK','','',1.0,1.0,1);
    %S232('dropall');
elseif (ephoneflag & ~nLocs)        % Use earphone filters only
    S232('PD1clrsched',1);
    S232('PD1nstrms',1,2,0);
    S232('PD1resetDSP',1,hex2dec('FFF'));
    S232('dropall');
    %Make connections for left ear
    S232('PD1addsimp',1,S232('DSPout',0),S232('DAC',0));        %DSPout to DAC0
    S232('PD1specIB',1,S232('IB',0),S232('DSPin',0));           %IB to DSPin
    %Make connections for right ear
    S232('PD1addsimp',1,S232('DSPout',1),S232('DAC',1));
    S232('PD1specIB',1,S232('IB',1),S232('DSPin',1));
    %Load left      
    S232('pushf',ephonefilt_left,length(ephonefilt_left));
    S232('PreLoadRaw',1,S232('DSPid',0),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    S232('dropall');
    %Load right
    S232('pushf',ephonefilt_right,length(ephonefilt_right));
    S232('PreLoadRaw',1,S232('DSPid',1),'MONO','STACK','','',TDT.ephonescale,1.0,1);
    S232('dropall');
elseif ~ephoneflag
    S232('PD1clrsched',1);
    S232('PD1nstrms',1,2,0);
    S232('PD1specIB',1,S232('IB',0),S232('DAC',0));
    S232('PD1specIB',1,S232('IB',1),S232('DAC',1));
end

% mute attenuators
S232('PA4mute',1);
S232('PA4mute',2);

%Set MII parameters
mii_us_per_sample = 10;             % microsecond per sample
mii_separation = 100;               % only take events separated by 100 samples (i.e., 1 ms)


ISI = XStimParams.test_ISI;
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000));         % correct for ITD filtlength
if(exist1('H.fig') & get(H.ephoneuseit,'Value'))    % correct for ephone filtlength
    ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));
end

%Loop to make the stimuli we'll use
if nLocs
    remtrials = nLocs * nABLs;
else
    remtrials = nILDs * nABLs * nITDs;
end

set(H.buildplay,'String','Building Stimuli');
set(H.remreps,'String',num2str(1));
set(H.remtrials,'string',num2str(remtrials));


%Add in the leading silent period
stim1 =  [zeros(1,silence_lead_pts) stim1(:)'];
stim2 = [zeros(1,silence_lead_pts) stim2(:)'];
%Add in the trailing silent period
stim1 =  [stim1 zeros(1,silence_trail_pts)];
stim2 = [stim2 zeros(1,silence_trail_pts)];

trialnum = 0;
set(H.buildplay,'BackgroundColor','red');
for iABL = 1:nABLs          % loop through ABLs
    for iITD = 1:nITDs          % loop through ITDs
        for iILD = 1:nILDs      % loop through ILDs
            %Check for pause by user
            if pause_check  return; end
            
            %%%%%%%%% apply ILD (half to each side)
            ILD = XStimParams.ilds(iILD);
            trial_left  = stim1 * (10^(-ILD/(2*20)));
            trial_right = stim2 * (10^(ILD/(2*20)));
            
            %%%%%%%%% Apply ITD filtering
            itdleft = 0; itdright = 0;
            ITD = round(XStimParams.itds(iITD));
            if(ITD < 0)
                itdright = abs(XStimParams.itds(iITD));
            elseif(ITD > 0)
                itdleft = abs(XStimParams.itds(iITD));
            end
            eval(['load ' FN.itd_path 'itdfilt' num2str(itdleft)]);
            eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
            itd_filt_left = itd_filt_left/max(abs(itd_filt_left));
            eval(['load ' FN.itd_path 'itdfilt' num2str(itdright)]);
            eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
            itd_filt_right = itd_filt_right/max(abs(itd_filt_left));
            trial_left = conv(trial_left,itd_filt_left);
            trial_right = conv(trial_right,itd_filt_right);
            
            if ephoneflag
                trial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
                trial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
            end
            
            % scaleFactor
            trial_left = trial_left * TDT.scaleFactor;
            trial_right = trial_right * TDT.scaleFactor;
            
            if nLocs            % loop through Locations
                trial_left = [trial_left zeros(1,TDT.hrtffiltlen)];
                trial_right = [trial_right zeros(1,TDT.hrtffiltlen)];
                for iLoc = 1:nLocs      % loop through Locations
                    
                    temp = conv(trial_left, hrTF1(iLoc,:));
                    % save to disk
                    trialnum=trialnum+1;
                    locs(:,trialnum) = XStimParams.locations(:,iLoc);
                    ilds(trialnum) = ILD;
                    itds(trialnum) = ITD;
                    abls(trialnum) = ABL;
                    S232('push16',temp,length(temp));
                    S232('qpop16',BUF.L1);
                    fname = ['LRstim.left_' num2str(trialnum)];
                    evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
                            [' ''' FN.temp_stim_path fname ''' '] ',0);'];
                    eval(evalstr);
                    temp = conv(trial_right, hrTF2(iLoc,:));
                    S232('push16',temp,length(temp));
                    S232('qpop16',BUF.R1);
                    fname = ['LRstim.right_' num2str(trialnum)];
                    evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
                            [' ''' FN.temp_stim_path fname ''' '] ',0);'];
                    eval(evalstr);
                end
            else               %%%%%% save stimuli to disk
                trialnum=trialnum+1;
                ilds(trialnum) = ILD;
                itds(trialnum) = ITD;
                locs = [];
                S232('push16',trial_left,length(trial_left));
                S232('qpop16',BUF.L1);
                fname = ['LRstim.left_' num2str(trialnum)];
                evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
                        [' ''' FN.temp_stim_path fname ''' '] ',0);'];
                eval(evalstr);
                S232('push16',trial_right,length(trial_right));
                S232('qpop16',BUF.R1);
                fname = ['LRstim.right_' num2str(trialnum)];
                evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
                        [' ''' FN.temp_stim_path fname ''' '] ',0);'];
                eval(evalstr);
            end         % iLoc
            abls(trialnum) = abs(XStimParams.abls(iABL));
        end             % iILD
    end                 % iITD
end                 % iABL
nTrials = trialnum;

%%%% finished making stimuli
% play out one sound
%Set up MII
m100x( C_.INIT );
m110dx( C_.INIT );
m110dx( C_.CLOCK, mii_us_per_sample);
m110dx( C_.MODE, M_.PST );
m110dx( C_.START);
m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse

nReps = XStimParams.numreps;
set(H.buildplay,'BackgroundColor','yellow');

%Put up raster plot
position = [900 50 300 900];
xvals = 1;
yvals = 1:nTrials;
xtext = '';
ytext = 'stim #';
xoffset = 1;
yoffset = 1;
tot_dur = XStimParams.silence_lead + DUR + 50;
stim_dur = XStimParams.silence_lead + DUR;
startcount = XStimParams.silence_lead;
endcount = XStimParams.silence_lead + DUR;
row_space = 5;
col_space = .2*tot_dur;

[H.rasterfig,H.rastaxes] = makerasterfig(position,...
    xvals, yvals,...
    xtext, ytext,...
    xoffset, yoffset,...
    nReps,...
    tot_dur, stim_dur,...
    startcount, endcount,...
    row_space, col_space);

% for plotting
finalspikematrix = zeros(1,nTrials);

%Begin loop over reps and trials   
set(H.buildplay,'String','Playing Stimuli');
set(H.buildplay,'BackgroundColor','yellow');
set(H.remreps,'String',num2str(nReps));
repnum = 1;
datamatrix = [];

if(exist1('H.LRsoundsfig') & get(H.recorddata,'Value'))
    update_dataFN;
end
while (repnum <= nReps)
    %Randomize the stimuli
    randseq = randperm(nTrials);   
    trialnum = 1;
    spikes_trial = [];
    repnum_trial = [];
    Nspikes = [];
    itds_trial = [];
    ilds_trial = [];
    abls_trial = [];
    el_trial = [];
    az_trial = [];
    iTrial_trial =[];
    while (trialnum <= nTrials+1)
        %Check for pause by user
        if pause_check  return; end
        
        if(trialnum <= nTrials)
            iTrial = randseq(trialnum);         % the randomly picked FN# to be loaded
            fname = ['LRstim.left_' num2str(iTrial)];
            evalstr = ['S232(''disk2dama16'',BUF.L1,' ...
                    [' ''' stimuli_dir fname ''' '] ',0);'];
            eval(evalstr);
            fname = ['LRstim.right_' num2str(iTrial)];
            evalstr = ['S232(''disk2dama16'',BUF.R1,' ...
                    [' ''' stimuli_dir fname ''' '] ',0);'];
            eval(evalstr);
        end
        
        %Wait till PD1 is finished
        pause(DUR/1000+.1);     % to bypass PD1status bug
        S232('PD1stop',1);
        %Check for pause by user
        if pause_check  return; end
        % set attenuators
        S232('PA4atten',1,abls(iTrial)-20);
        S232('PA4atten',2,abls(iTrial)-20);
        
        if trialnum > 1
            iTrial = randseq(trialnum-1);
            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, 1000);                    %Take 1000 spikes max
            spikes = spikes(find(spikes ~= 0));                 %Get clock events that are spikes
            ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
            if ~isempty(ind)
                spikes = spikes(ind);
                spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
                itds_trial = [itds_trial;itds(iTrial) * ones(size(spikes))];
                ilds_trial = [ilds_trial;ilds(iTrial) * ones(size(spikes))];
                abls_trial = [abls_trial;abls(iTrial) * ones(size(spikes))];
                if nLocs
                    el_trial = [el_trial; locs(1,trialnum)* ones(size(spikes))];
                    az_trial = [az_trial; locs(2,trialnum)* ones(size(spikes))];
                end
                repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
                Nspikes = [Nspikes; ones(length(spikes),1)*length(spikes)];
                iTrial_trial = [iTrial_trial; iTrial * ones(size(spikes))];
            end
        end
        
        if trialnum <= nTrials
            S232('seqplay',BUF.playspecbuf);
            S232('PD1arm',1);
            %Set up MII
            m100x( C_.INIT );
            m110dx( C_.INIT );
            m110dx( C_.CLOCK, mii_us_per_sample);
            m110dx( C_.MODE, M_.PST );
            pause(ISI/1000);
            if (trialnum <= nTrials)        %Start clock & Send pulse: PD1 GO!
                m110dx( C_.START);
                m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
            end
        end
        
        if(trialnum > 1)
            ind0 = find(spikes < 100*XStimParams.silence_lead & spikes > 0);            % spikes preceding sound onset
            ind1 = find(spikes > 100*XStimParams.silence_lead & spikes <= (XStimParams.silence_lead + DUR)*100);   % spikes during sound
            evoked =  length(ind1) - length(ind0);
            finalspikematrix(iTrial) = finalspikematrix(iTrial) + evoked;
        end
        
        remtrials = nTrials - trialnum;
        set(H.remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end %end loop over trials
    %Plot all the spikes from trials
    plotraster(H.rasterfig,...
        spikes_trial,...
        ones(size(spikes_trial)),...
        iTrial_trial,...
        repnum,...
        tot_dur, stim_dur,...
        nReps,...
        xoffset, yoffset,...
        row_space, col_space);
    
    %Record Data
    if get(H.recorddata,'Value')
        if nLocs
            datamatrix = [datamatrix; [Nspikes spikes_trial repnum_trial itds_trial ilds_trial abls_trial el_trial az_trial]];
        else
            datamatrix = [datamatrix; [Nspikes spikes_trial repnum_trial itds_trial ilds_trial abls_trial]];
        end
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
    plot(1:nTrials,finalspikematrix/repnum,'b-*','LineWidth',1.5);
    set(H.axes,'Color','black');
    xlabel('Trial #'); ylabel('Spike Rate (spikes/stim)');
    title('LR test');
    %%%%%
    
    set(H.remreps,'String',num2str(nReps - repnum));
    repnum = repnum + 1;
    pause(0);
end %end loop over reps


set(H.buildplay,'String','Build/Play status');
set(H.exitLRsounds,'Visible','on');
set(H.resetLRsounds,'Visible','on');

% increment testnum
if(exist1('H.LRsoundsfig') & get(H.recorddata,'Value'))
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
while (exist1('H.LRsoundsfig') & get(H.pauseLRsounds,'Value'))
    pause(0);
    if(~exist1('H.LRsoundsfig')) return; end         
    set(H.exitLRsounds,'Visible','on');
    set(H.resetLRsounds,'Visible','on');
    if(exist1('H.LRsoundsfig') & get(H.resetLRsounds,'Value') == 1)
        set(H.resetLRsounds,'Value',0);
        set(H.pauseLRsounds,'Value',0);
        Reset_ITD_decorr;  flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end