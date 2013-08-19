function [] = Engage_AdaptedSpace()

%% DECLARE GLOBALS
global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_AdaptedSpace

%*******************************************************************************
%*******************************************************************************
 colors = [ ...
         179 199 255; ...
        200 150 255; ...
        0    0  255; ...
        216 41  0; ...
        255 199 179;...
        255 150 200;...
        255   0    0;...
        199 255 179;...
        200 255 150;...
        0   255 0]/255;

fclose all;

%% SCALE FACTOR SETTINGS (High Dynamic Range?)
if XStimParams.HiDynamicRange
    scaleFactor = 60;
else
    scaleFactor = TDT.scaleFactor;
end

% reset stim_type to file (#9)
set(H.stim_type,'Value',9);

XStimParams.curr_ITD = str2num(get(H.AdaptedSpace_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.AdaptedSpace_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.AdaptedSpace_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.AdaptedSpace_ISI,'String'));
XStimParams.numreps = str2num(get(H.AdaptedSpace_numreps,'String'));
XStimParams.reset_flag = 0;

% load filter coeffs & other params for ABLalone test
if XStimParams.ABLalone_flag
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

%Specify DAMA buffers
clear BUF
BUF.AdaptL				= 1; % adaptor (left ch)
BUF.ProbeL              = 2; % probe (left ch)
BUF.AdaptR              = 3; % adaptor (right ch)
BUF.ProbeR				= 4; % probe (right ch)
BUF.playseq_L1          = 5;
BUF.playseq_R1          = 6;
BUF.playspec1           = 7;

%Make play sequence buffers
S232('allot16',BUF.playseq_L1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.AdaptL);
S232('make',1,BUF.ProbeL);
S232('make',2,1);
S232('make',3,0);
S232('qpop16',BUF.playseq_L1);

S232('allot16',BUF.playseq_R1,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.AdaptR);
S232('make',1,BUF.ProbeR);
S232('make',2,1);
S232('make',3,0);
S232('qpop16',BUF.playseq_R1);

%Make play specification buffer
S232('allot16',BUF.playspec1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L1);
S232('make',1,BUF.playseq_R1);
S232('make',2,0);
S232('qpop16',BUF.playspec1);

% DETERMINE TOTAL TRIAL DURATION (ADAPTOR+PROBE) ... FIX THIS
TOTAL_DUR = ADAPTOR_DUR + PROBE_DUR;

% Fully-Cued Test
Npts_adaptor = (ADAPTOR_DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen + TDT.ephonefiltlen;
Npts_probe = (PROBE_DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen + TDT.ephonefiltlen;
Npts_totalplay = Npts_adaptor + Npts_probe;
XStimParams.test_type = 'AdaptedSpace FC'; str1 = 'fc';
disp('This is a FULLY CUED AdaptedSpace test')

S232('allot16',BUF.AdaptL,Npts_adaptor);
S232('allot16',BUF.ProbeL,Npts_probe);
S232('allot16',BUF.AdaptR,Npts_adaptor);
S232('allot16',BUF_ProbeR,Npts_probe);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,Npts_totalplay);

%Get Earphone filters
if FN.HRTFfiletype(6) == 1
    ephonefname = [FN.ephone_path FN.ephone2];
    ephonefilt_left  = (mtlrch(ephonefname,1))';
    ephonefilt_right = (mtlrch(ephonefname,2))';
else
    eval(['load -mat ' FN.ephone_path FN.ephone2]);
    ephonefilt_left  = TF1;
    ephonefilt_right = TF2;
    clear TF1 TF2 dir
end

%Load Earphone filters
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
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000)); 				%correct for ITD filtlength
ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));			%correct for ephonefilt

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
if get(H.AdaptedSpace_recorddata,'Value')
    tempstr = ['    ' str1 ' AdaptedSpace-test: ' num2str(size(XStimParams.locations,2)) ' locations;     atten: ' num2str(abs(ABL))];
    update_diary
end

Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_AdaptedSpace_' str1 ' Temp_params;'])
clear Temp_params str1

if isempty(FN.space_std) | ~isempty(strfind(FN.space_std,'eq'))        % FN not yet picked
    [FN.space_std,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.std preferred)');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',[FN.space_path FN.space_std]);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
end
%%%%%%%
if FN.HRTFfiletype(1,2) == 1
    hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space_std]));
elseif FN.HRTFfiletype(1,2) == 2
    dir = 0;
    eval(['load -mat ' FN.space_path FN.space_std]);
    TF1_space = TF1;    TF2_space = TF2;
    hrtfdirmat = dir;
    clear dir TF1 TF2
else
    disp(['space HRTFfiletype incorrect'])
    return
end

clear locind
for locnum = 1:size(XStimParams.locations,2)
    locind(locnum) = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
end

% update display
setInfo_AdaptedSpace;

%%%%% Loop to make the adaptor stimuli
remreps = 1;
set(H.AdaptedSpace_status,'String','Status: building stimuli');
set(H.AdaptedSpace_status,'BackgroundColor','yellow');
set(H.AdaptedSpace_status,'ForegroundColor','blue');
set(H.AdaptedSpace_remreps,'String',num2str(remreps));
repnum = 1;
nLocs = size(XStimParams.locations,2);
finalspikematrix = zeros(1,numtrials);

% clear out temporary directory
delete( [FN.temp_stim_path '\*.*'])

%% MAKE ADAPTOR

% Check for pause by user
if pause_check    return; end

% make adaptor (by type) here FIX THIS!

adapt_left = adapt_left(:)' - mom(adapt_left,1);

% trailing pad stim length of short stims
if length(adapt_left)/TDT.Fs *1000 < DUR
    npts = ceil(DUR*30 - length(adapt_left));
    adapt_left = [adapt_left zeros(1,npts)];
end

% normalize to ACPower
if ~XStimParams.HiDynamicRange
    adapt_left = adapt_left / mom(adapt_left,2);
end
adapt_right = adapt_left;
adapt_left0 = adapt_left;
adapt_right0 = adapt_right;

for iLoc = 1:nLocs        %Apply HRTF filtering
    if pause_check    return; end        %Check for pause by user

    if FN.HRTFfiletype(1,2) == 1
        hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locind(iLoc))-1);
        hrtf_right = mtlrch([FN.space_path FN.space_std],2*locind(iLoc));
    else
        hrtf_left = TF1_space(locind(iLoc),:);
        hrtf_right = TF2_space(locind(iLoc),:);
    end
    str1 = 'FC';
    
    adapt_left = conv(adapt_left0, hrtf_left);
    adapt_right = conv(adapt_right0, hrtf_right);
    
    % remove DC offset
    adapt_left = adapt_left - mean(adapt_left);
    adapt_right = adapt_right - mean(adapt_right);
    
    % adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
    % presented through DSPs (*.eq) at 0,0 (broadband)
    adapt_left = adapt_left * TDT.hrtf_Lfactor;
    adapt_right = adapt_right * TDT.hrtf_Rfactor;
    
    % scale stimuli 3/7/07
    % note scale factor changes with XStimParams.HiDynamicRange
    adapt_left = adapt_left * scaleFactor;
    adapt_right = adapt_right * scaleFactor;
    
    % pad with zeros for ephonefilters
    adapt_left = [adapt_left zeros(1,TDT.ephonefiltlen)];
    adapt_right = [adapt_right zeros(1,TDT.ephonefiltlen)];
    
    % check that adapt_left and adapt_right have Npts_adaptor???? FIX THIS
    
    % save stimuli to disk
    S232('push16',adapt_left,length(adapt_left));
    S232('qpop16',BUF.AdaptL);
    fname = ['ADAPTOR_' str1 '.left_' ...
        num2str(hrtfdirmat(1,locind(iLoc))) ...
        '_' num2str(hrtfdirmat(2,locind(iLoc)))];
    evalstr = ['S232(''dama2disk16'',BUF.AdaptL,' [' ''' FN.temp_stim_path fname ''' '] ',0);'];
    eval(evalstr);
    
    S232('push16',adapt_right,length(adapt_right));
    S232('qpop16',BUF.AdaptR);
    fname = ['ADAPTOR_' str1 '.right_' ...
        num2str(hrtfdirmat(1,locind(iLoc))) ...
        '_' num2str(hrtfdirmat(2,locind(iLoc)))];
    evalstr = ['S232(''dama2disk16'',BUF.AdaptR,' [' ''' FN.temp_stim_path fname ''' '] ',0);'];
    eval(evalstr);
    
    % set up for next trial
    remstims = nLocs - iLocs + 1;
    set(H.AdaptedSpace_remstims,'String',num2str(remstims));
    pause(0);
end     % iLoc

%% Make Test Probe

% Check for pause by user
if pause_check    return; end

% make probe (by type) here FIX THIS!

probe_left = probe_left(:)' - mom(probe_left,1);

% trailing pad stim length of short stims
if length(probe_left)/TDT.Fs *1000 < DUR
    npts = ceil(DUR*30 - length(probe_left));
    probe_left = [probe_left zeros(1,npts)];
end

% normalize to ACPower
if ~XStimParams.HiDynamicRange
    probe_left = probe_left / mom(probe_left,2);
end
probe_right = probe_left;
probe_left0 = probe_left;
probe_right0 = probe_right;

if FN.HRTFfiletype(1,2) == 1
    hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locind(iLoc))-1);
    hrtf_right = mtlrch([FN.space_path FN.space_std],2*locind(iLoc));
else
    hrtf_left = TF1_space(locind(iLoc),:);
    hrtf_right = TF2_space(locind(iLoc),:);
end
str1 = 'FC';

probe_left = conv(probe_left0, hrtf_left);
probe_right = conv(probe_right0, hrtf_right);

% remove DC offset
probe_left = probe_left - mean(probe_left);
probe_right = probe_right - mean(probe_right);

% adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
% presented through DSPs (*.eq) at 0,0 (broadband)
probe_left = probe_left * TDT.hrtf_Lfactor;
probe_right = probe_right * TDT.hrtf_Rfactor;

% scale stimuli 3/7/07
% note scale factor changes with XStimParams.HiDynamicRange
probe_left = probe_left * scaleFactor;
probe_right = probe_right * scaleFactor;

% pad with zeros for ephonefilters
probe_left = [probe_left zeros(1,TDT.ephonefiltlen)];
probe_right = [probe_right zeros(1,TDT.ephonefiltlen)];

% check that probe_left and probe_right have Npts_adaptor???? FIX THIS

% save stimuli to disk
S232('push16',probe_left,length(probe_left));
S232('qpop16',BUF.ProbeL);
fname = ['PROBE_' str1 '.left_' ...
    num2str(hrtfdirmat(1,locind(iLoc))) ...
    '_' num2str(hrtfdirmat(2,locind(iLoc)))];
evalstr = ['S232(''dama2disk16'',BUF.ProbeL,' [' ''' FN.temp_stim_path fname ''' '] ',0);'];
eval(evalstr);

S232('push16',probe_right,length(probe_right));
S232('qpop16',BUF.ProbeR);
fname = ['PROBE_' str1 '.right_' ...
    num2str(hrtfdirmat(1,locind(iLoc))) ...
    '_' num2str(hrtfdirmat(2,locind(iLoc)))];
evalstr = ['S232(''dama2disk16'',BUF.ProbeR,' [' ''' FN.temp_stim_path fname ''' '] ',0);'];
eval(evalstr);

set(H.AdaptedSpace_remstims,'String',num2str(0));
pause(0);
%%%%%%%%%%%%%%%%%%% finished making stimuli

%Begin playing sounds   
set(H.AdaptedSpace_status,'String','Status: playing stimuli');
set(H.AdaptedSpace_status,'BackgroundColor','green');
set(H.AdaptedSpace_remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.AdaptedSpacefig') & get(H.AdaptedSpace_recorddata,'Value'))
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
%loop for reps
    trialnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    FN_trial = [];
    Nspikes = [];
while (exist1('H.AdaptedSpacefig') & (repnum <= XStimParams.numreps))
    tic
    
    % create a randomized sequence of conditions (Adaptor Loc x Probe Scale)
    % top row = adaptor locations
    % bottom row = probe spls
    tmp1=repmat(1:nLocs,1,nProbeScales); tmp2=[];
    for i=1:nProbeScales
        tmp2=[tmp2 i*ones(1,nLocs)];
    end
    tmpSEQ=[tmp1; tmp2];
    nSEQ=nLocs*nProbeScales;
    SEQ=zeros(2,nSEQ);
    randSEQ=randperm(nSEQ);
    for j=1:nSEQ
        SEQ(:,j)=tmpSEQ(:,randSEQ(j));
    end
    
    % loop through sequence of conditions
    for iSEQ = 1:nSEQ
        
        %Check for pause by user
        if pause_check    return; end
        
        %%% load ADAPTOR stimuli
        
        % get this trials adaptor loc
        iLoc=SEQ(1,iSEQ);
        
        %% push left file onto stack & check filesize
        fname = ['ADAPTOR_' str1 '.left_' ...
            num2str(hrtfdirmat(1,locind(iLoc))) ...
            '_' num2str(hrtfdirmat(2,locind(iLoc)))];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
        eval(evalstr);
        % check if file is correct length (FIX THIS)
        
        %% pop to DAMA
        eval(['S232(''qpop16'''  ',BUF.AdaptL);']);
        
        %% push right file onto stack & check filesize
        fname = ['ADAPTOR_' str1 '.right_' ...
            num2str(hrtfdirmat(1,locind(iLoc))) ...
            '_' num2str(hrtfdirmat(2,locind(iLoc)))];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
        eval(evalstr);
        % check if file is correct length (FIX THIS)
        
        %% pop to DAMA
        eval(['S232(''qpop16'''  ',BUF.AdaptR);']);
        
        %% LOAD PROBE STIMULI !!!
        %(FIX THIS)
        trial_scale=ProbeScales(SEQ(2,iSEQ));
        
        %% push left file onto stack & check filesize
        fname = ['PROBE_' str1 '.left_' ...
            num2str(hrtfdirmat(1,locind(ProbeLocIdx))) ...
            '_' num2str(hrtfdirmat(2,locind(ProbeLocIdx)))];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
        eval(evalstr);
        % check if file is correct length (FIX THIS)
        
        %% pop to DAMA
        eval(['S232(''qpop16'''  ',BUF.ProbeL);']);
        
        %% push right file onto stack & check filesize
        fname = ['PROBE_' str1 '.right_' ...
            num2str(hrtfdirmat(1,locind(ProbeLocIdx))) ...
            '_' num2str(hrtfdirmat(2,locind(ProbeLocIdx)))];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
        eval(evalstr);
        % check if file is correct length (FIX THIS)
        
        %% pop to DAMA
        eval(['S232(''qpop16'''  ',BUF.ProbeR);']);
        
        S232('seqplay',BUF.playspec1);
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
        
        while(S232('PD1status',1)) usec_delay(1000);  end
        S232('PD1stop',1);
        
        %Stop the m110 and get spiketimes
        m110dx( C_.STOP);
        spikes = m110dx( C_.DATA, round(XStimParams.curr_stimdur*2)); 			% Take 2*XStimParams.curr_stimdur spikes max
        ind = find(spikes ~= 0); 						% Get clock events that are spikes
        spikes = spikes(ind);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
        if(exist1('H.AdaptedSpacefig') & ~isempty(spikes))
            spikes = spikes(ind);
            spikes_trial = [spikes_trial; spikes/(1000/mii_us_per_sample)];
            AdaptEL_trial = [AdaptEL_trial; hrtfdirmat(1,locind(iLoc))* ones(size(spikes))];
            AdaptAZ_trial = [AdaptAZ_trial; hrtfdirmat(2,locind(iLoc))* ones(size(spikes))];
            repnum_trial = [repnum_trial; repnum * ones(size(spikes))];
            ProbeScale_trial = [ProbeScale_trial; trial_scale * ones(size(spikes))];F
            Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
        end
        ind = (iLoc-1)*nFiles + iFile;
        finalspikematrix(ind) = finalspikematrix(ind) + length(spikes);
        if pause_check    return; end
        
        remtrials = numtrials - trialnum;
        set(H.AdaptedSpace_remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end             % iLocs
    
    %Record Data
    if(exist1('H.AdaptedSpacefig') & get(H.AdaptedSpace_recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial AdaptEL_trial AdaptAZ_trial ProbeScale_trial]];
        record_data3(XStimParams,datamatrix);
    end
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.AdaptedSpacefig') & ~exist1('H.AdaptedSpace_finalspikeratefig'))
        H.AdaptedSpace_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','AdaptedSpace Test Spike Rate Plot',...
            'NumberTitle','off');
        H.AdaptedSpace_spikeaxes = axes;
    end
    figure(H.AdaptedSpace_finalspikeratefig); hold off
    %plotdiam1(XStimParams.locations, interimspikerate);
    for iLoc = 1:nLocs
       ind = [1:nFiles] + (iLoc-1)*nFiles;
       plot(1:nFiles,finalspikematrix(ind),'color',colors(rem(iLoc,10)+1,:))
       hold on
    end
    set(H.AdaptedSpace_spikeaxes,'Color','black');
    xlabel('file#'); ylabel('spikerate'); title(['Rep # ' num2str(repnum)]);
            
    remreps = XStimParams.numreps - repnum;
    set(H.AdaptedSpace_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.AdaptedSpace_finalspikeratefig); hold off
set(H.AdaptedSpace_finalspikeratefig,'Name','Final Plot for AdaptedSpace Test');
    for iLoc = 1:nLocs
       ind = [1:nFiles] + (iLoc-1)*nFiles;
       plot(1:nFiles,finalspikematrix(ind),'color',colors(rem(iLoc,10)+1,:))
       hold on
    end
set(H.AdaptedSpace_spikeaxes,'Color','black');
xlabel('file#'); ylabel('spikerate');
title('good stuff - finished!', 'FontSize',8);

set(H.AdaptedSpace_status,'String','Status: results');
set(H.AdaptedSpace_status,'BackgroundColor','blue');

set(H.exitAdaptedSpace,'Visible','on');
set(H.resetAdaptedSpace,'Visible','on');

% increment test number
if(exist1('H.AdaptedSpacefig') & get(H.AdaptedSpace_recorddata,'Value'))
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
while (exist1('H.AdaptedSpacefig') & get(H.pauseAdaptedSpace,'Value'))
    pause(0);
    if(~exist1('H.AdaptedSpacefig')) return; end         
    set(H.exitAdaptedSpace,'Visible','on');
    set(H.resetAdaptedSpace,'Visible','on');
    if(exist1('H.AdaptedSpacefig') & get(H.resetAdaptedSpace,'Value') == 1)
        set(H.resetAdaptedSpace,'Value',0);
        set(H.pauseAdaptedSpace,'Value',0);
        Reset_AdaptedSpace;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_AdaptedSpace;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end