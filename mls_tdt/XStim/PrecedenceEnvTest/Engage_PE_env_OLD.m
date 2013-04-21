function [] = Engage_PE_env()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI
global FigParams

% Engage_PE_env

%*******************************************************************************
%	The Precedence Effect AM-Test operation
%*******************************************************************************
% uses *.eq HRTFfiles
%       
% Jun08_09 changed default scaling to .5, so that total level will
% approximate one source

rand('state',sum(100*clock));
rampPts = TDT.Fs *2.5/1000;             % stim made separately to 2.5 ms ramps
stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);

% scaleFactor   %%%%%%%debug   *10 to use as spike trigger
scaleFactor = TDT.scaleFactor;  

set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

% calc stim duration
DUR = XStimParams.curr_stimdur;
nPts = DUR*round(TDT.Fs/1000);
nEpochmsecs = XStimParams.epoch_duration(1);       % duration of each epoch at one SPL
nEpochPts = round(nEpochmsecs * TDT.Fs/1000);
nEpochs = round(nPts / nEpochPts);
nPts = nEpochs * nEpochPts;
DUR = nPts / round(TDT.Fs/1000);
XStimParams.curr_stimdur = DUR;

if XStimParams.space_flag
    while exist1([FN.space_path FN.space_eq]) ~=2
        [FN.space_eq,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select Fully-cued HRTF File *.eq preferred');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space_eq]);
        end
        FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);
    end
    disp('This is a FULLY-CUED PE AMstimuli test')
else
    while exist1([FN.space_path FN.ablequal_eq]) ~=2
        [FN.ablequal_eq,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select ABLequal FC HRTF File  *.eq preferred');
        if(FN.space_path ~= 0)
            set(H.ABLequalfile,'String',[FN.space_path FN.ablequal_eq]);
        end
        FN.HRTFfiletype(7,1) = testHRTFfiletype(FN.space_path, FN.ABLequal_eq);
    end
    disp('This is an ABLequal PE AMstimuli test')
end    

%Put parameters into XStimParams
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
%XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.reset_flag = 0;

%Specify DAMA buffers
clear BUF
BUF.snd1			= 1;
BUF.snd2			= 2;
BUF.zeros			= 3;
BUF.playseq0        = 4;        % BEWARE: don't change this without changing below
BUF.playseq1		= 5;        % BEWARE: don't change this without changing below
BUF.playseq2		= 6;        % BEWARE: don't change this without changing below
BUF.playspec1		= 7;

%Make play sequence buffers
S232('allot16',BUF.playseq0,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.zeros);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq0);

S232('allot16',BUF.playseq1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.snd1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq1);

S232('allot16',BUF.playseq2,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.snd2);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq2);

%Make play specification buffer
S232('allot16',BUF.playspec1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq1);
S232('make',1,BUF.playseq2);
S232('make',2,0);
S232('qpop16',BUF.playspec1);

%Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
silence_len = (XStimParams.silence_lead * round(TDT.Fs/1000));
%Add a piece of silence after stimulus 
silence_len2 = (XStimParams.silence_trail * round(TDT.Fs/1000));

%Make Stimulus buffers
nPts_totalplay = (silence_len2 + silence_len + nPts) + 2* TDT.hrtffiltlen;
S232('dpush',nPts_totalplay);
S232('value',0);
S232('allot16',BUF.zeros,nPts_totalplay);
S232('qpop16',BUF.zeros);

S232('dpush',nPts_totalplay);
S232('value',0);
S232('allot16',BUF.snd1,nPts_totalplay);
S232('qpop16',BUF.snd1);

S232('dpush',nPts_totalplay);
S232('value',0);
S232('allot16',BUF.snd2,nPts_totalplay);
S232('qpop16',BUF.snd2);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,nPts_totalplay);

S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('FFF'));           % was '0xFFF'
S232('dropall');

%Make connections for left ear
sources(1) = S232('DSPout',0);          % snd1
sources(2) = S232('DSPout',2);          % snd2
scales(1) = 1;
scales(2) = 1;
S232('PD1addmult',1,sources,scales,2,S232('DAC',0));

%Make connections for right ear
sources(1) = S232('DSPout',1);      % snd1
sources(2) = S232('DSPout',3);      % snd2
scales(1) = 1;
scales(2) = 1;
S232('PD1addmult',1,sources,scales,2,S232('DAC',1));

S232('PD1addsimp',1,S232('IREG',0),S232('DSPin',0));
S232('PD1addsimp',1,S232('IREG',0),S232('DSPin',1));
S232('PD1specIB',1,S232('IB',0),S232('IREG',0));

S232('PD1addsimp',1,S232('IREG',1),S232('DSPin',2));
S232('PD1addsimp',1,S232('IREG',1),S232('DSPin',3));
S232('PD1specIB',1,S232('IB',1),S232('IREG',1));

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4mute',1);	
S232('PA4mute',2);	

%ISI = XStimParams.test_ISI;

%Get HRTF index
XStimParams.locations = GUI.locations1';
if size(XStimParams.locations,2) ~= 1
    set(H.pickerfig,'Color', [.1 .2 .8]);
    set(H.picker_error,'visible','on');
    disp('PE AM test requires one and only one picked location')
    pause;
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end

%%%% load HRTFs
dir = 0;
if XStimParams.space_flag
    tempFN = FN.space_eq;
    filetype = FN.HRTFfiletype(1,1);
else
    tempFN = FN.ablequal_eq;
    filetype = FN.HRTFfiletype(7,1);
end
if filetype == 2
    eval(['load -mat ' FN.space_path tempFN]);
elseif filetype == 1
    dir = sph2dbl(mtlrdir([FN.space_path tempFN]));
end
XStimParams.HRTF_FN = tempFN;

% find locations for 2 sources NOTE:: OFFSETS here are for loc2 vs loc1
% (loc2 = loc1 + offset)
% (not, as in Two Source, either source vs mid-array)
offset_el = XStimParams.offset_el(2);
offset_az = XStimParams.offset_az(2);

loc1 = max(find(dir(1,:) == XStimParams.locations(1) & dir(2,:) == XStimParams.locations(2)));
loc2 = max(find(dir(1,:) == (XStimParams.locations(1) + offset_el) & dir(2,:) == (XStimParams.locations(2)) + offset_az));
loc2use = [loc1 loc2];
str1 = ['  ' num2str(XStimParams.locations(1)) ' el, ' num2str(XStimParams.locations(2)) ' az, vs: ' ...
        num2str(XStimParams.locations(1) + offset_el) ' el, ' num2str(XStimParams.locations(2) + offset_az) ' az'];

if isempty(loc1)  || isempty(loc2)
    error('offset not possible, repick')
    return;
end

%%%%%%% debug - set loc2 HRTFs ==0
%TF1(loc2,:) = TF1(loc2,:)*0;
%TF2(loc2,:) = TF2(loc2,:)*0;

%%%% new
%load DSPs for sound#1 HRTF filtering
%Load left   
S232('dropall');
S232('pushf',TF1(loc1,:),length(TF1(loc1,:)));
S232('PreLoadRaw',1,S232('DSPid',0),'MONO','STACK','','',TDT.ephonescale,1.0,1);
%Load right
S232('pushf',TF2(loc1,:),length(TF2(loc1,:)));
S232('PreLoadRaw',1,S232('DSPid',1),'MONO','STACK','','',TDT.ephonescale,1.0,1);

%load DSPs for sound#2 HRTF filtering
%Load left   
S232('dropall');
S232('pushf',TF1(loc2,:),length(TF1(loc2,:)));
S232('PreLoadRaw',1,S232('DSPid',2),'MONO','STACK','','',TDT.ephonescale,1.0,1);
%Load right
S232('pushf',TF2(loc2,:),length(TF2(loc2,:)));
S232('PreLoadRaw',1,S232('DSPid',3),'MONO','STACK','','',TDT.ephonescale,1.0,1);
%%%%%%

% write diary
if get(H.recorddata,'Value')
    tempstr = ['     PE_env   atten: ' num2str(abs(XStimParams.curr_ABL)) 'dB'];
    update_diary
end
clear hrtf_*

S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);
pause(.1);

% increment testnumber
if(exist1('H.PE_envfig') & get(H.recorddata,'Value'))
    update_dataFN;
end

% Begin playing sounds
tic
maxSpikes = ceil(nPts_totalplay/10);    % max # spikes to download from M110d
tempPts = rampPts:nPts-rampPts;
cmap = hsv(length(XStimParams.cF)*2);
set(H.status,'String','Status: Playing Stimuli');
set(H.status,'BackgroundColor','green');
set(H.remreps,'String',num2str(XStimParams.numreps));
irep = 1;
datamatrix = [];

% how many tests??
numtests=0;
n_cF = length(XStimParams.cF);
n_delay = length(XStimParams.delay);
for iLoc=1:2
    for icF = 1:n_cF
        for idelay = 1:n_delay
            numtests= numtests+1;
        end
    end
end
% decorr env, decorr carrier
% identical env, decorr carrier
% identical env, identical carrier one source in RF
% identical env, identical carrier one source in anti-RF
for icF = 1:n_cF
    numtests= numtests+4;
end

usedBank = zeros(XStimParams.numreps, numtests);        % keep track of rows of randBank used by each test
if exist('e:\kip\MATLAB\mls_tdt\XStim\PrecedenceEnvTest\randBank.mat','file')
    load 'e:\kip\MATLAB\mls_tdt\XStim\PrecedenceEnvTest\randBank.mat'
else
    randBank = ceil(rand(XStimParams.numreps,5) * 10000000);
end

%loop for reps (of blocks)
%Htemp = figure;
while (exist1('H.PE_envfig') & (irep <= XStimParams.numreps))
    % make params for this rep
    [AMparams{irep} usedBank] = makeAMparams(randBank,usedBank);
    %figure(Htemp);
    %pcolor(usedBank); caxis([0 2]);
    testOrder = randperm(numtests);
    set(H.remreps,'String',num2str(XStimParams.numreps - irep));

    % loop for tests (= 1 block)
    spike_raster_matrix = cell(numtests,1);
    itest = 1;
    while itest <= numtests
        
        % make & load next stim & load DSPs
        if irep <= XStimParams.numreps
            set(H.remtrials,'String',num2str(numtests - itest));
            [Noise1, Noise2] = filteredstimuliAM_01(AMparams{irep}(testOrder(itest),:), 1);
            % define locs 1 and 2 as locs for snd 1 & 2, respectively
            %loc1 = loc2use(AMparams{irep}(testOrder(itest),11));
            %if AMparams{irep}(testOrder(itest),12)
            %    loc2 = loc2use(AMparams{irep}(testOrder(itest),12));
            %else
            %    loc2 = 0;
            %end
            % scale to AC power == 1 * scaleFactor
            %if AMparams{irep}(testOrder(itest),11) == AMparams{irep}(testOrder(itest),12)    % single source
            %    temp = std(Noise1(tempPts) + Noise2(tempPts));
            %    Noise1 = Noise1/temp * scaleFactor;
            %    Noise2 = Noise2/temp * scaleFactor;
            %else
            Noise1 = Noise1/std(Noise1(tempPts)) * scaleFactor;
            Noise2 = Noise2/std(Noise2(tempPts)) * scaleFactor;
                %end            
            
            
            %%%%%%%% debug   set Noise2 ==0
            %Noise1 = Noise1*0;
          
            S232('dropall');
            S232('push16',[zeros(1,silence_len) Noise1 zeros(1,silence_len2+ 255*2)],nPts_totalplay);
            S232('qpop16',BUF.snd1);
            S232('push16',[zeros(1,silence_len) Noise2 zeros(1,silence_len2+ 255*2)],nPts_totalplay);
            S232('qpop16',BUF.snd2);
            
            if 0  %%% old way
                %load DSPs for sound#1 HRTF filtering
                %Load left   
                S232('dropall');
                S232('pushf',TF1(loc1,:),length(TF1(loc1,:)));
                S232('PreLoadRaw',1,S232('DSPid',0),'MONO','STACK','','',TDT.ephonescale,1.0,1);
                %Load right
                S232('pushf',TF2(loc1,:),length(TF2(loc1,:)));
                S232('PreLoadRaw',1,S232('DSPid',1),'MONO','STACK','','',TDT.ephonescale,1.0,1);
                
                %load DSPs for sound#2 HRTF filtering
                %Load left   
                S232('dropall');
                S232('pushf',TF1(loc2,:),length(TF1(loc2,:)));
                S232('PreLoadRaw',1,S232('DSPid',2),'MONO','STACK','','',TDT.ephonescale,1.0,1);
                %Load right
                S232('pushf',TF2(loc2,:),length(TF2(loc2,:)));
                S232('PreLoadRaw',1,S232('DSPid',3),'MONO','STACK','','',TDT.ephonescale,1.0,1);
            else      %%%%%%%% new
                %Make play specification buffer
                
                if AMparams{irep}(testOrder(itest),12) == 0 & AMparams{irep}(testOrder(itest),11) == 2
                    S232('dropall');
                    S232('dpush',1);
                    S232('make',0,AMparams{irep}(testOrder(itest),11)+4);     %%% BEWARE +4 makes these correspond to assignments near top               
                    S232('qpoppart16',BUF.playspec1,1);
                    
                    S232('dpush',1);
                    S232('make',0,AMparams{irep}(testOrder(itest),12)+4);     %%% BEWARE +4 makes these correspond to assignments near top
                    S232('qpoppart16',BUF.playspec1,0);
                else                
                    S232('dropall');
                    S232('dpush',1);
                    S232('make',0,AMparams{irep}(testOrder(itest),11)+4);     %%% BEWARE +4 makes these correspond to assignments near top               
                    S232('qpoppart16',BUF.playspec1,0);                 % RF
                    
                    S232('dpush',1);
                    S232('make',0,AMparams{irep}(testOrder(itest),12)+4);     %%% BEWARE +4 makes these correspond to assignments near top
                    S232('qpoppart16',BUF.playspec1,1);                 % antiRF
                end
            end
            
            %%%%%% debug
            disp(['loc1:   ' num2str(AMparams{irep}(testOrder(itest),11)) '    loc#2:   ' num2str(AMparams{irep}(testOrder(itest),12))])
            %pause(.1)
            
        end             % end of making next stim & loading DSPs
        
        while toc < (XStimParams.curr_stimdur + 500)/1000     end
        tic

        S232('seqplay',BUF.playspec1);
        S232('PD1arm',1);
        %Set up MII & send trigger
        m100x( C_.INIT );
        m110dx( C_.INIT );
        m110dx( C_.CLOCK, mii_us_per_sample);
        m110dx( C_.MODE, M_.PST );
        
        m110dx( C_.START);
        m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
        
        %Wait for TOC
        while toc <= (XStimParams.curr_stimdur+300)/1000     end
 
        S232('PD1stop',1);
        %Stop the m110 and get spikes
        m110dx( C_.STOP);
        spikes = m110dx( C_.DATA, maxSpikes);
        %ind = find(spikes ~= 0); 						% Get clock events that are spikes
        spikes = spikes(spikes ~= 0);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
        if ~isempty(spikes) 
            spikes = spikes(ind);
            spikes_trial = spikes/(1000/mii_us_per_sample);
            EL_trial = dir(1,loc2use(AMparams{irep}(testOrder(itest),11)))* ones(size(spikes));
            AZ_trial = dir(2,loc2use(AMparams{irep}(testOrder(itest),11)))* ones(size(spikes));
            irep_trial = irep * ones(size(spikes));
            itest_trial = testOrder(itest) * ones(size(spikes));
            Nspikes = length(spikes) * ones(size(spikes));
            spike_raster_matrix{testOrder(itest)} = spikes_trial;
            
            %Record Data
            if get(H.recorddata,'Value')
                datamatrix = [datamatrix;[Nspikes spikes_trial irep_trial itest_trial EL_trial AZ_trial]];
                record_data3(XStimParams,datamatrix, AMparams);
            end
        end
        
        if pause_check  return; end
        pause(0)
        itest = itest+1;
    end         % itest (= 1 block)
    
    %%%%%%%%%%%
    %Plot Spike Raster
    if ~exist1('H.finalspikerasterfig')
        H.finalspikerasterfig = figure('Position',[700 20 550 500],...
            'Name','PE_env Test Spike Raster',...
            'NumberTitle','off');
        H.spikeaxes = axes;
        set(H.spikeaxes,'Color','black');
        hold on
        axis([0 XStimParams.curr_stimdur + XStimParams.silence_lead + XStimParams.silence_trail 0 XStimParams.numreps*(1+numtests)])
        set(H.spikeaxes,'Ytick',0:XStimParams.numreps:XStimParams.numreps*numtests);
        set(H.spikeaxes,'YtickLabel',0:1:numtests);
        
        ylabel('test#')
        xlabel('time (ms)')
        title('spike raster plot')
        
        plot([1 1]*XStimParams.silence_lead,[0 XStimParams.numreps*(1+numtests)],'g')
        plot([1 1]*(XStimParams.silence_lead+XStimParams.curr_stimdur),[0 XStimParams.numreps*(1+numtests)],'r')
    end
    figure(H.finalspikerasterfig);
    for itest = 1:numtests
        icF = find(XStimParams.cF == AMparams{1}(itest,7));
        temp = spike_raster_matrix{itest}(spike_raster_matrix{itest}<XStimParams.silence_lead | spike_raster_matrix{itest}>XStimParams.silence_lead + +XStimParams.curr_stimdur);
        plot(temp,ones(size(temp))*itest*XStimParams.numreps-irep+1,'.')
        temp = spike_raster_matrix{itest}(spike_raster_matrix{itest}>XStimParams.silence_lead & spike_raster_matrix{itest}<XStimParams.silence_lead + +XStimParams.curr_stimdur);
        plot(temp,ones(size(temp))*itest*XStimParams.numreps-irep+1,'.','color', cmap(icF+1,:))
    end
        
    %Plot spike envelope plot
    Plot_PE_env(datamatrix,AMparams);
    
    irep = irep + 1;
    if pause_check  return; end
    pause(0);
end 									%end loop over reps

S232('PA4mute',1);	
S232('PA4mute',2);	

set(H.status,'String','Status: Results');
set(H.exitPE_env,'Visible','on');
set(H.resetPE_env,'Visible','on');

% increment test number
if(exist1('H.PE_envfig') & get(H.recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
while (exist1('H.PE_envfig') & get(H.pausePE_env,'Value'))
    pause(0);
    if(~exist1('H.PE_envfig')) return; end         
    set(H.exitPE_env,'Visible','on');
    set(H.resetPE_env,'Visible','on');
    if(exist1('H.PE_envfig') & get(H.resetPE_env,'Value') == 1)
        set(H.resetPE_env,'Value',0);
        set(H.pausePE_env,'Value',0);
        Reset_PE_env;   flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_PE_env;   flag=1;
        return;
    end
end

if XStimParams.reset_flag
    flag = 1;
    XStimParams.reset_flag = 0;
end