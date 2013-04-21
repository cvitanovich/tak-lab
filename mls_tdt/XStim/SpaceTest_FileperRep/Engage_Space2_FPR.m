function [] = Engage_space2_FPR()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

% Engage_space2_FPR
%*******************************************************************************
%	The Space Test operation 3/15/07
% eliminated double buffering
% does not write to files
% maximum duration about 30 secs
% FPR = filePerRep  (differing files for each "rep" - played through nLocs
% and over nReps (true reps of each FN)
% set up to read *.noi files (floats)
%*******************************************************************************
colors = [ ...
        .7 .7 1; ...
        .8 .6 1; ...
        0  0  1; ...
        .8 .2 0; ...
        1 .8 .7;...
        1 .6 .8;...
        1 0  0;...
        .8 1 .7;...
        .8 1 .6;...
        0 1 0];

XStimParams.test_type = 'space_FPR FC';        str1 = 'FC';
disp('This is a FULLY CUED SPACE test')

fclose all;
eval(['delete ' FN.temp_stim_path '*.*;']);

if XStimParams.HiDynamicRange
    scaleFactor = 60;
else
    scaleFactor = TDT.scaleFactor;
end

% reset stim_type to file (#9)
set(H.stim_type,'Value',9);

pts_per_msec = round(TDT.Fs/1000);

set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

XStimParams.curr_ITD = str2num(get(H.space_FPR_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.space_FPR_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.space_FPR_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space_FPR_ISI,'String'));
XStimParams.numreps = str2num(get(H.space_FPR_numreps,'String'));
XStimParams.reset_flag = 0;

%Specify DAMA buffers
clear BUF
BUF.L1  		    = 1;        % for playing sounds
BUF.R1  		    = 2;        % for playing sounds
BUF.silence_lead    = 3;
BUF.silence_trail   = 4;
BUF.playseq_L		= 5;
BUF.playseq_R		= 6;
BUF.playspec		= 7;

S232('PD1stop',1);
S232('PD1clear',1);
S232('dropall');
S232('trash')

%Make play sequence buffers
S232('allot16',BUF.playseq_L,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.silence_lead);
S232('make',1,1);
S232('make',2,BUF.L1);
S232('make',3,1);
S232('make',4,BUF.silence_trail);
S232('make',5,1);
S232('make',6,0);
S232('qpop16',BUF.playseq_L);

S232('allot16',BUF.playseq_R,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.silence_lead);
S232('make',1,1);
S232('make',2,BUF.R1);
S232('make',3,1);
S232('make',4,BUF.silence_trail);
S232('make',5,1);
S232('make',6,0);
S232('qpop16',BUF.playseq_R);

%Make play specification buffer
S232('allot16',BUF.playspec,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L);
S232('make',1,BUF.playseq_R);
S232('make',2,0);
S232('qpop16',BUF.playspec);

% check to be sure stimFNs chosen
if isempty(FN.FPR)         %Stimulus from file
    set(H.space_FPR_initFN,'value',1)
    setinfo_space_FPR
end
nFiles = length(FN.FPR);
for iFile=1:nFiles
    D = dir([FN.stim_path2 FN.FPR{iFile}]);
    B(iFile) = D.bytes/4;
end
DUR = max1(B) / pts_per_msec;
XStimParams.curr_stimdur = DUR;

% make leading and lagging silence buffers
S232('allot16',BUF.silence_lead,(XStimParams.silence_lead * pts_per_msec));
S232('dpush',(XStimParams.silence_lead * pts_per_msec));
S232('value',0);
S232('qpop16',BUF.silence_lead);

S232('allot16',BUF.silence_trail,(XStimParams.silence_trail * pts_per_msec)+TDT.hrtffiltlen);
S232('dpush',(XStimParams.silence_trail * pts_per_msec)+TDT.hrtffiltlen);
S232('value',0);
S232('qpop16',BUF.silence_trail);

% Make Stimulus buffers
S232('allot16',BUF.L1,DUR * pts_per_msec);
S232('allot16',BUF.R1,DUR * pts_per_msec);

% set up PD1
S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
Npts_totalplay = (DUR + XStimParams.silence_lead + XStimParams.silence_trail) * pts_per_msec + TDT.hrtffiltlen;
S232('PD1npts',1,Npts_totalplay);

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

%Set MII parameters
mii_us_per_sample = 10; 			%microsecond per sample
mii_separation = 100; 				%only take events separated by 100 samples (i.e., 1 ms)

ITD = XStimParams.curr_ITD;
if(abs(ITD) > 250) return; end

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4mute',1);
S232('PA4mute',2);

ISI = XStimParams.test_ISI;

% save XStimParams for this test
if get(H.space_FPR_recorddata,'Value')
    tempstr = [str1 ' space_FPR-test: ' num2str(length(XStimParams.locations)) ' locations;     atten: ' num2str(abs(ABL))];
    update_diary
end
Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_space_FPR_' str1 ' Temp_params;'])
clear Temp_params str1

% check whether hrtf Files picked
if XStimParams.space_flag
    if isempty(FN.space_eq)       % FN not yet picked
        [FN.space_eq,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.eq preferred)');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space_eq]);
        end
        set(H.spacefilepb,'Value',0);
    end
    t1 = 1;   tempPATH = FN.space_path;   tempFN = FN.space_eq;    str1 = 'FC';
end   
FN.HRTFfiletype(t1,1) = testHRTFfiletype(tempPATH, tempFN);

%%%%%%%
if FN.HRTFfiletype(t1,1) == 1
    hrtfdirmat = sph2dbl(mtlrdir([tempPATH tempFN]));
    for iLoc = 1:size(hrtfdirmat,2)
        hrTF1(iLoc,:) = mtlrch([tempPATH tempFN],(2*iLoc-1));
        hrTF2(iLoc,:) = mtlrch([tempPATH tempFN],(2*iLoc));
    end
elseif FN.HRTFfiletype(t1,1) == 2
    dir = 0;
    eval(['load -mat ' tempPATH tempFN]);
    hrtfdirmat = dir; clear dir
    hrTF1 = TF1;    hrTF2 = TF2; clear TF*
else
    disp(['HRTFfiletype incorrect'])
    return
end

%Get HRTF indices
XStimParams.locations = GUI.locations1';
if ~length(XStimParams.locations)
    set(H.pickerfig,'Color', [.1 .2 .8]);
    set(H.picker_error,'visible','on');
    pause;
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end
nLocs = size(XStimParams.locations,2);
clear locind
for iLoc = 1:nLocs
    locind(iLoc) = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,iLoc) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,iLoc)));
end

% update display
setInfo_space_FPR;

%%%%% play out one sound to get rid of duration-change artifact
set(H.space_FPR_status,'BackgroundColor','yellow');
set(H.space_FPR_status,'ForegroundColor','blue');
set(H.space_FPR_status,'String','Status: loading stimulus');

%%%% load file
if strcmp('mat',FN.FPR{1}(end-2:end))
    eval(['load ' FN.stim_path2 FN.FPR{1}])
    trial_left  = stim;
else
    fid = fopen([FN.stim_path2 FN.FPR{1}],'r');
    trial_left = fread(fid,inf,'float');
    fclose(fid);
end
trial_left = trial_left(:)' - mom(trial_left,1);
trial_left = trial_left / mom(trial_left,2);

%scale stimuli 3/7/07
trial_left = trial_left * scaleFactor;

% trailing pad stim length of short stims
if length(trial_left)/TDT.Fs *1000 < DUR
    npts = ceil(DUR*pts_per_msec - length(trial_left));
    trial_left = [trial_left zeros(1,npts)];
end

% load to buffers
S232('push16',trial_left,DUR * pts_per_msec);
S232('qpop16',BUF.L1);
S232('push16',trial_left,DUR * pts_per_msec);
S232('qpop16',BUF.R1);

%Begin playing sounds   
set(H.space_FPR_status,'String','Status: playing silent PRE_stimulus');
set(H.space_FPR_status,'BackgroundColor','green');
set(H.space_FPR_status,'ForegroundColor','white');

% increment testnumber
if(exist1('H.space_FPRfig') & get(H.space_FPR_recorddata,'Value'))
    update_dataFN;
end

S232('seqplay',BUF.playspec);
S232('PD1arm',1);

%Set up MII
m100x( C_.INIT );
m110dx( C_.INIT );
m110dx( C_.CLOCK, mii_us_per_sample);
m110dx( C_.MODE, M_.PST );
%Start clock
m110dx( C_.START);
tic
m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse

% Wait for these buffers to finish playing
t = toc;
while t < (Npts_totalplay/TDT.Fs)+ISI/1000
    t = toc;
end

m110dx( C_.STOP);
spikes = m110dx( C_.DATA, 100); 			% Take 100 spikes max

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

%loop for reps
datamatrix = [];
repnum = 1;
numtrials = nLocs * nFiles;
finalspikematrix = zeros(1,numtrials);

set(H.space_FPR_status,'String','Status: playing stimuli');
set(H.space_FPR_status,'BackgroundColor','green');
set(H.space_FPR_status,'ForegroundColor','white');

% set up raster figure
flag_Raster = 1;
if (strcmp(FN.FPR{1}(1:4),'Cond'))      % conditioner_Probe
    FPRtype = 1;
    yvals = fliplr([2.5 5 10 20 40 80 160 320 640 1280 2560]);
    ytext = 'interval (ms)';
elseif (strcmp(FN.FPR{1}(1:2),'ND'))    % Nagel_Doupe
    FPRtype = 2;
    yvals = fliplr(1:nFiles);
    ytext = 'FileNumber';
else
    flag_Raster = 0;
end
if flag_Raster
    position = [900 50 300 900];
    xvals = 1;
    xtext = '';
    xoffset = 1;
    yoffset = 1;
    tot_dur = XStimParams.silence_lead + DUR + 50;
    stim_dur = XStimParams.silence_lead + DUR;
    startcount = XStimParams.silence_lead;
    endcount = XStimParams.silence_lead + DUR;
    row_space = 5;
    col_space = .2*tot_dur;
    
    [H.rasterfig,H.rastaxes] = makerasterfig(position,...
        xvals,...
        yvals,...
        xtext,...
        ytext,...
        xoffset,...
        yoffset,...
        XStimParams.numreps,...
        tot_dur,...
        stim_dur,...
        startcount,...
        endcount,...
        row_space,...
        col_space);
end

%Randomize stimFNs
FN.FPR = FN.FPR(randperm(nFiles));

while (exist1('H.space_FPRfig') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli
    trialnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    FN_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    for iFile = 1:nFiles
        randseq = randperm(nLocs);
        set(H.space_FPR_status,'BackgroundColor','yellow');
        set(H.space_FPR_status,'ForegroundColor','blue');
        set(H.space_FPR_status,'String','Status: loading stimulus');
        
        %%%% load file
        if strcmp('mat',FN.FPR{iFile}(end-2:end))
            eval(['load ' FN.stim_path2 FN.FPR{iFile}])
            trial_left  = stim;
        else
            fid = fopen([FN.stim_path2 FN.FPR{iFile}],'r');
            trial_left = fread(fid,inf,'float');
            fclose(fid);
        end
        trial_left = trial_left(:)' - mom(trial_left,1);
        trial_left = trial_left / mom(trial_left,2);
        
        % trailing pad stim length of short stims
        if length(trial_left)/TDT.Fs *1000 < DUR
            npts = ceil(DUR*pts_per_msec - length(trial_left));
            trial_left = [trial_left zeros(1,npts)];
        end

        %scale stimuli 3/7/07
        trial_left = trial_left * scaleFactor;
                
        % load to buffers
        S232('push16',trial_left,DUR * pts_per_msec);
        S232('qpop16',BUF.L1);
        S232('push16',trial_left,DUR * pts_per_msec);
        S232('qpop16',BUF.R1);
        
        set(H.space_FPR_status,'String','Status: playing stimuli');
        set(H.space_FPR_status,'BackgroundColor','green');
        set(H.space_FPR_status,'ForegroundColor','white');
        
        for iLoc = 1:nLocs
            %Check for pause by user
            if pause_check    return; end
            
            % load HRTFs
            %Load left      
            S232('pushf',hrTF1(locind(randseq(iLoc)),:),255);
            S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
            %Load right
            S232('pushf',hrTF2(locind(randseq(iLoc)),:),255);
            S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.hrtf_scale,1.0,1);
            
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
            m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
            tic
            
            pause(DUR/1000+.4);
            
            %while(S232('PD1status',1)) usec_delay(1000);  end
            %S232('PD1stop',1);
                        
            %Stop the m110 and get spikes
            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, XStimParams.curr_stimdur * 2); 			    % Take (2* dur in msec) spikes max
            ind = find(spikes ~= 0); 						% Get clock events that are spikes
            spikes = spikes(ind);
            ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
            if(exist1('H.space_FPRfig') & ~isempty(spikes)) 
                spikes = spikes(ind);
                spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
                EL_trial = [EL_trial;hrtfdirmat(1,locind(randseq(iLoc)))* ones(size(spikes))];
                AZ_trial = [AZ_trial;hrtfdirmat(2,locind(randseq(iLoc)))* ones(size(spikes))];
                FN_trial = [FN_trial; iFile* ones(size(spikes))];
                repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
                Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
            end
            ind = (iLoc-1)*nFiles + iFile;
            finalspikematrix(ind) = finalspikematrix(ind) + length(spikes);
            
            
            % raster
            if(exist1('H.rasterfig') & ~isempty(spikes)) & flag_Raster
                switch FPRtype
                    case 1      % Conditioner_Probe
                        ind = strfind(FN.FPR{iFile},'_');
                        if strfind(FN.FPR{iFile}(ind(1)+1:ind(2)-1),'p')
                            FNnum = 11;
                        else
                        FNnum = find(yvals == str2num(FN.FPR{iFile}(ind(1)+1:ind(2)-1)));
                    end
                    case 2      % Nagel_Doupe
                        ind = strfind(FN.FPR{iFile},'rep');
                        %FNnum = find(yvals == str2num(FN.FPR{iFile}(ind(1)+3:end)));
                        FNnum = iFile;
                end
                plotraster(H.rasterfig,...
                    spikes/(1000/mii_us_per_sample),...
                    ones(size(spikes)),...
                    ones(size(spikes))*FNnum,...
                    repnum,...
                    tot_dur,...
                    stim_dur,...
                    XStimParams.numreps,...
                    xoffset,...
                    yoffset,...
                    row_space,...
                    col_space);
                drawnow
            end
            
            %Record Data
            if(exist1('H.space_FPRfig') & get(H.space_FPR_recorddata,'Value'))
               record_data3(XStimParams,[Nspikes spikes_trial repnum_trial EL_trial AZ_trial FN_trial]);
            end
            
            if pause_check    return; end
            
            remtrials = numtrials - trialnum;
            set(H.space_FPR_remtrials,'String',num2str(remtrials));
            trialnum = trialnum + 1;
        end             % iLocs
    end                 % iFiles
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.space_FPRfig') & ~exist1('H.space_FPR_finalspikeratefig'))
        H.space_FPR_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','Space Test Spike Rate Plot',...
            'NumberTitle','off');
        H.space_FPR_spikeaxes = axes;
    end
    figure(H.space_FPR_finalspikeratefig)
    for iLoc = 1:nLocs
        ind = [1:nFiles] + (iLoc-1)*nFiles;
        if nLocs ==1
            plot(1:nFiles,interimspikerate(ind),'co')
        else
            plot(1:nFiles,interimspikerate(ind),'o','color',colors(rem(iLoc,10)+1,:))
        end
        hold on
    end
    set(H.space_FPR_spikeaxes,'Color','black');
    xlabel('file#'); ylabel('spikerate'); title(['Rep # ' num2str(repnum)]);
    
    %Record Data
    if(exist1('H.space_FPRfig') & get(H.space_FPR_recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial FN_trial]];
        record_data3(XStimParams,datamatrix);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.space_FPR_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
    
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.space_FPR_finalspikeratefig)
set(H.space_FPR_finalspikeratefig,'Name','Final Plot for Space Test');
for iLoc = 1:nLocs
    ind = [1:nFiles] + (iLoc-1)*nFiles;
    if nLocs == 1
        plot(1:nFiles,finalspikematrix(ind),'yo')
    else
        plot(1:nFiles,finalspikematrix(ind),'o','color',colors(rem(iLoc,10)+1,:))
    end
    hold on
end
set(H.space_FPR_spikeaxes,'Color','black');
xlabel('file#'); ylabel('spikerate');
title('good stuff - finished!', 'FontSize',8);

set(H.space_FPR_status,'String','Status: results');
set(H.space_FPR_status,'BackgroundColor','blue');

set(H.exitspace_FPR,'Visible','on');
set(H.resetspace_FPR,'Visible','on');

% increment test number
if(exist1('H.space_FPRfig') & get(H.space_FPR_recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end


%%%%%%%%%%%
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
while (exist1('H.space_FPRfig') & get(H.pausespace_FPR,'Value'))
    pause(0);
    if(~exist1('H.space_FPRfig')) return; end         
    set(H.exitspace_FPR,'Visible','on');
    set(H.resetspace_FPR,'Visible','on');
    if(exist1('H.space_FPRfig') & get(H.resetspace_FPR,'Value') == 1)
        set(H.resetspace_FPR,'Value',0);
        set(H.pausespace_FPR,'Value',0);
        Reset_space;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_space_FPR;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end