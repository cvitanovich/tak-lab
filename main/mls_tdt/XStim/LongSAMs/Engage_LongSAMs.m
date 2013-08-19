function [] = Engage_LongSAMs(tempStruct)

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_LongSAMs

%*******************************************************************************
%	The LongSAMs Test operation
% should use HRIRs with earphone equalization - usually called *.eq
%
%       play2_record2b_SAMB.dll or
%       play2_record2b_SAMc measured on 1/28/09 introduces a 97.9 msec
%       time advance into test. That is, the spike-clock starts 97.9 msec after
%       the sound - so the spiketimes appear earlier than they really are


%       remeasured the time problem on 8/16/10, by using the envelope click
%       commented out below
%           found a delay of 1.83 to 1.96 msec (that is stim starts ca 2
%           msec after spike collection)
%
%

%       - the MII clock is slow: 1 sec on the TDT == 998.348 msec on the
%       m110d (Jan 08, 2009)

%       altered Oct 20, 2009 to allow double modulation, where
%       mod_freq(2) is the (main) slower rate
%       currently, both modulations are weighted equally
%*******************************************************************************
XStimParams.test_type = 'LongSAMs FC';     str1 = 'FC';
disp('This is a FULLY CUED LongSAMs test')

stimuli_dir = FN.temp_stim_path;
fclose all;

if 0
if XStimParams.HiDynamicRange
    scaleFactor = 60;
else
    scaleFactor = TDT.scaleFactor;
end
end

% reset stim_type to BroadBand
set(H.stim_type,'Value',8);

XStimParams.curr_ABL = round(str2num(get(H.LongSAMs_ABL,'String')));
set(H.LongSAMs_ABL,'String',num2str(XStimParams.curr_ABL));
ABL = abs(XStimParams.curr_ABL);
XStimParams.curr_stimdur = str2num(get(H.LongSAMs_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.LongSAMs_ISI,'String'));
XStimParams.numreps = str2num(get(H.LongSAMs_numreps,'String'));
XStimParams.reset_flag = 0;

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

%Get the HRTF indices
XStimParams.locations = GUI.locations1';
if ~length(XStimParams.locations)
    set(H.pickerfig,'Color', [.1 .2 .8]);
    set(H.picker_error,'visible','on');
    pause;
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end

% save XStimParams for this test
if get(H.LongSAMs_recorddata,'Value')
    if strcmp(XStimParams.mod_type, 'dBTone')
        str2 = 'dB'
    else
        str2 = '';
    end
    tempstr = ['    ' str1 str2 ' LongSAMs-test: ' num2str(XStimParams.mod_freq(1)) ' Hz @' num2str(100*XStimParams.mod_depth(1))  '% x ' num2str(tempStruct.nCycles_perRep) ' cycles x ' num2str(XStimParams.numreps) ' reps    atten: ' num2str(abs(ABL))];
    update_diary
end

Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_LongSAMs_' str1 ' Temp_params;'])
clear Temp_params str1

if isempty(FN.space_eq) | ~isempty(strfind(FN.space_eq,'std'))        % FN not yet picked
    [FN.space_eq,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.eq preferred)');
    if(FN.space_path ~= 0)
        set(H.LongSAMsfile,'String',[FN.space_path FN.space_eq]);
    end
    set(H.LongSAMsfilepb,'Value',0);
    FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);
end
%%%%%%%
if FN.HRTFfiletype(1,1) == 1
    hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space_eq]));
elseif FN.HRTFfiletype(1,1) == 2
    dir = 0;
    eval(['load -mat ' FN.space_path FN.space_eq]);
    TF1_LongSAMs = TF1;    TF2_LongSAMs = TF2;
    hrtfdirmat = dir;
    clear dir TF1 TF2
else
    disp(['LongSAMs HRTFfiletype incorrect'])
    return
end
clear locind
for locnum = 1:size(XStimParams.locations,2)
    locind(locnum) = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
end

% update display
setInfo_LongSAMs;

% update tempStruct for play2_record2b_SAMB
tempStruct.mdepth = XStimParams.mod_depth(1);
tempStruct.latten = abs(XStimParams.curr_ABL)-20;
tempStruct.ratten = abs(XStimParams.curr_ABL)-20;
tempStruct.SamplingPeriod = 1000000/TDT.Fs;
tempStruct.buf_pts = TDT.Fs / XStimParams.mod_freq(2);
tempStruct.nCycles_perRep = str2num(get(H.LongSAMstestmod_nCycles_perRep,'string'));
tempStruct.nptsTotalPlay = tempStruct.nCycles_perRep * tempStruct.buf_pts;
tempStruct.ISI = XStimParams.test_ISI;
tempStruct.record_spikes = 1;

% make envelope(s)
tt = [0:1/tempStruct.buf_pts:1-1/tempStruct.buf_pts];
Env = sin(2*pi*tt-pi/2);
switch XStimParams.LongSAMs_flag
    case 1
    case 2
        tempFreq = XStimParams.mod_freq(1) / XStimParams.mod_freq(2);
        Env1 = sin(2*pi*tt*tempFreq-pi/2);
        Env = (Env + Env1)/2;
    case 3                                  % if env1 is 1/16 Hz
        Env2 = sin(2*pi*tt*2-pi/2 +1);      % 1/8 Hz
        Env3 = sin(2*pi*tt*4-pi/2 +2);      % 1/4 Hz
        Env4 = sin(2*pi*tt*8-pi/2 +3);      % 1/2 Hz
        Env = (Env + Env2 + Env3 + Env4)/4;
    case 4                                  % if env1 is 1/16 Hz
        Env2 = sin(2*pi*tt*2-pi/2 +1);      % 1/8 Hz
        Env3 = sin(2*pi*tt*4-pi/2 +2);      % 1/4 Hz
        Env4 = sin(2*pi*tt*8-pi/2 +3);      % 1/2 Hz
        Env5 = sin(2*pi*tt*16-pi/2 +2.5);   % 1 Hz
        Env6 = sin(2*pi*tt*32-pi/2 +3.5);   % 2 Hz
        Env7 = sin(2*pi*tt*64-pi/2 +4);     % 4 Hz
        Env = (Env + Env2 + Env3 + Env4 + Env5 + Env6 + Env7)/7;
end

if strcmp(XStimParams.mod_type, 'dBTone')
    Env = (Env +1) /2;
    Env = 10 .^((XStimParams.mod_dB(1)*Env-XStimParams.mod_dB(1))/20);
else
    Env = Env * (XStimParams.mod_depth(1)/2) + (1 - XStimParams.mod_depth(1)/2);
end


% temporary for testing
%Env =  zeros(size(Env));
%Env(1:10) = ones(1,10)*2;

Env = single(Env);

%Begin playing sounds   
set(H.LongSAMs_status,'String','Status: playing stimuli');
set(H.LongSAMs_status,'BackgroundColor','green');
set(H.LongSAMs_remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
buffcycle = 1;
datamatrix = [];
numlocs = length(locind);
finalspikematrix = zeros(1,numlocs);

% increment testnumber
if(exist1('H.LongSAMsfig') & get(H.LongSAMs_recorddata,'Value'))
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
%loop for reps
while (exist1('H.LongSAMsfig') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli
    randseq = randperm(numlocs);   
    locnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    tic
    while (exist1('H.LongSAMsfig') & (locnum <= numlocs) )
        %Check for pause by user
        if pause_check    return; end
        
        % load stimuli
        HRTF_L = single(TF1_LongSAMs(locind(randseq(locnum)),:));
        HRTF_R = single(TF2_LongSAMs(locind(randseq(locnum)),:));
        
        %Set up MII (done in .dll)
        %m100x( C_.INIT );
        %m110dx( C_.INIT );
        %m110dx( C_.CLOCK, mii_us_per_sample);
        %m110dx( C_.MODE, M_.PST );
        
        while toc < XStimParams.test_ISI/1000 
            if pause_check    return; end
        end

        %Start clock
        tic
        %m110dx( C_.START);

        play2_record2b_SAMc(tempStruct, HRTF_L(:)', HRTF_R(:)', Env);
        %play2_record2b_SAMb(tempStruct, HRTF_L(:)', HRTF_R(:)');
        while toc < (XStimParams.curr_stimdur + 1000)/1000 
            if pause_check    return; end
        end
        tic
        
        %Stop the m110 and get spikes
        m110dx( C_.STOP);
        nspikes = min1([XStimParams.curr_stimdur*2 50000]);
        spikes = m110dx( C_.DATA, nspikes); 			% Take nspikes max
        ind = find(spikes ~= 0); 						% Get clock events that are spikes
        spikes = spikes(ind);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
        if(exist1('H.LongSAMsfig') & ~isempty(spikes)) 
            spikes = spikes(ind);
            spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
            EL_trial = [EL_trial;hrtfdirmat(1,locind(randseq(locnum)))* ones(size(spikes))];
            AZ_trial = [AZ_trial;hrtfdirmat(2,locind(randseq(locnum)))* ones(size(spikes))];
            repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
            Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
        end
        finalspikematrix(randseq(locnum)) = finalspikematrix(randseq(locnum)) + length(spikes);
        if pause_check    return; end
    
        remtrials = numlocs - locnum;
        set(H.LongSAMs_remtrials,'String',num2str(remtrials));
        locnum = locnum + 1;
        pause(0);
    end %end loop over trials
    
    if 0
        %Plot Spike Rate Data
        interimspikerate = finalspikematrix/repnum;
        if(exist1('H.LongSAMsfig') & ~exist1('H.LongSAMs_finalspikeratefig'))
            H.LongSAMs_finalspikeratefig = figure('Position',[700 20 550 500],...
                'Name','LongSAMs Test Spike Rate Plot',...
                'NumberTitle','off');
            H.LongSAMs_spikeaxes = axes;
        end
        figure(H.LongSAMs_finalspikeratefig)
        plotdiam1(XStimParams.locations, interimspikerate);
        set(H.LongSAMs_spikeaxes,'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
        colorbar
    end
    
    %Record Data
    datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
    if(exist1('H.LongSAMsfig') & get(H.LongSAMs_recorddata,'Value'))
        DATA = record_data3(XStimParams,datamatrix);
    elseif ~get(H.LongSAMs_recorddata,'Value')
        DATA = record_data3(XStimParams,datamatrix,[],[],[],0);
    end
    
    % plot period histogram
    if(exist1('H.LongSAMsfig') & ~exist1('H.LongSAMs_PERIODfig'))
        H.LongSAMs_PERIODfig = figure('Position',[700 20 550 500],...
            'Name','LongSAMs Period Histogram',...
            'NumberTitle','off');
        H.LongSAMs_PERIODaxes = axes;
    end
    if length(DATA)>1
        pl_longSAMs(DATA,XStimParams.mod_freq(2),0);    % 3rd argin is 'overplot'
        if XStimParams.LongSAMs_flag == 2        % 2) double mod   1) single mod
            pl_longSAMs(DATA,XStimParams.mod_freq(1),1);
        end
    end
    remreps = XStimParams.numreps - repnum;
    set(H.LongSAMs_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

if 0
%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.LongSAMs_finalspikeratefig)
set(H.LongSAMs_finalspikeratefig,'Name','Final Plot for LongSAMs Test');
plotdiam1(XStimParams.locations, interimspikerate);
set(H.LongSAMs_spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
        ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))], 'FontSize',8);
colorbar
end
% update period histogram
figure(H.LongSAMs_PERIODfig)
title(['final period histogram: mod period = ' num2str(1000/ XStimParams.mod_freq(1)) ' msec  total Nspikes = ' num2str(sum(Nspikes))])

set(H.LongSAMs_status,'String','Status: results');
set(H.LongSAMs_status,'BackgroundColor','blue');

set(H.exitLongSAMs,'Visible','on');
set(H.resetLongSAMs,'Visible','on');

% increment test number
if(exist1('H.LongSAMsfig') & get(H.LongSAMs_recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
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
while (exist1('H.LongSAMsfig') & get(H.pauseLongSAMs,'Value'))
    pause(0);
    if(~exist1('H.LongSAMsfig')) return; end         
    set(H.exitLongSAMs,'Visible','on');
    set(H.resetLongSAMs,'Visible','on');
    if(exist1('H.LongSAMsfig') & get(H.resetLongSAMs,'Value') == 1)
        set(H.resetLongSAMs,'Value',0);
        set(H.pauseLongSAMs,'Value',0);
        Reset_LongSAMs;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_LongSAMs;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end