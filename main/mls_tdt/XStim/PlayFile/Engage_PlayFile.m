function [] = Engage_PlayFile()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_PlayFile

%*******************************************************************************
%	The PlayFile Test operation
% alterred version to eliminate double buffering
% should use HRIRs without earphone equalization - usually called *.std
%*******************************************************************************

stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);

if XStimParams.HiDynamicRange
    scaleFactor = 60;
else
    scaleFactor = TDT.scaleFactor;
end

if get(H.PlayFile_NoScalingpb,'value')
    scaleFlag = 0;
else
    scaleFlag = 1;
end

% reset non-allowed stim_type to BroadBand
if isempty(find([3 6 7 8 9] == get(H.stim_type,'Value')))
    set(H.stim_type,'Value',8);
    disp('Unsupported stimulus type for PlayFile test; reset to BROADBAND');
end
if get(H.stim_type,'Value') == 9
    if isempty(FN.stim)         %Stimulus from file
        [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select stimulus File');
    end
    D = dir([FN.stim_path FN.stim]);
    XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
end
set(H.PlayFile_DUR,'String',num2str(XStimParams.curr_stimdur));

XStimParams.curr_ITD = str2num(get(H.PlayFile_ITD,'String'));
XStimParams.curr_ILD = str2num(get(H.PlayFile_ILD,'String'));
XStimParams.curr_ABL = str2num(get(H.PlayFile_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.PlayFile_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.PlayFile_ISI,'String'));
XStimParams.numreps = str2num(get(H.PlayFile_numreps,'String'));
XStimParams.reset_flag = 0;
XStimParams.test_type = 'Play a File';     str1 = 'playFile';
disp('This is a FULLY CUED PlayFile test')

% save XStimParams for this test
Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_Space_' str1 ' Temp_params;'])
clear Temp_params

ITD = XStimParams.curr_ITD;
if(abs(ITD) > 250) return; end

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
ILD = XStimParams.curr_ILD;
ISI = XStimParams.test_ISI;
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000)); 				%correct for ITD filtlength
% update display
setInfo_PlayFile;

%Specify DAMA buffers
clear BUF
BUF.L1				= 1;
BUF.R1				= 2;
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

% Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
DUR_silence = XStimParams.silence_lead; 					%ms
silence_len = (DUR_silence * round(TDT.Fs/1000));
% Add a piece of silence after stimulus 
DUR_silence2 = XStimParams.silence_trail; 					%ms
silence_len2 = (DUR_silence2 * round(TDT.Fs/1000));

%Make Stimulus buffers
DUR = XStimParams.curr_stimdur;
Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000));

S232('allot16',BUF.L1,Npts_totalplay);
S232('allot16',BUF.R1,Npts_totalplay);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,Npts_totalplay);

% make schedule
S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1specIB',1,S232('IB',0),S232('DAC',0));
S232('PD1specIB',1,S232('IB',1),S232('DAC',1));

S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

%%%%% Loop to make the stimuli we'll use
remreps = 1;
set(H.PlayFile_status,'String','Status: building stimuli');
set(H.PlayFile_status,'BackgroundColor','blue');
set(H.PlayFile_remreps,'String',num2str(remreps));

finalspikematrix = zeros(1,1);
clear stim_from_file

set(H.PlayFile_status,'BackgroundColor','yellow');
%Check for pause by user
if pause_check    return; end

%Make the stimuli
switch get(H.stim_type,'Value')      
    case 3,
        %Use Concurrent ToneStack if specified: 10 frequencies
        tonestackfreqs = ...
            XStimParams.bandwidth(1):...
            round((XStimParams.bandwidth(2) - XStimParams.bandwidth(1))/9):...
            XStimParams.bandwidth(2);
        [trial_left,trial_right] = MakeToneStack(TDT.Fs,tonestackfreqs,XStimParams.curr_stimdur);
    case 6,
        %Narrowband Noise
        [trial_left] = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
        trial_right = trial_left;
    case 7,
        %Narrowband Noise Stack
        [trial_left] = MakeNBNoiseStack(TDT.Fs,XStimParams.curr_stimdur,XStimParams.spec_NBNoiseStackBands);
        trial_right = trial_left;
    case 8,
        %Broadband Noise
        [trial_left] = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        trial_right = trial_left;
    case 9,
        %Stimulus from file
        if(~exist1('stim_from_file'))
            fid = fopen([FN.stim_path FN.stim],'r');
            stim_from_file = fread(fid,inf,'float');
            fclose(fid);
        end
        [trial_left] = stim_from_file;
        if(size(trial_left,1) > 1)
            trial_left = trial_left';
        end
        trial_right = trial_left;
    otherwise,
        set(H.stim_type,'Value',8);
        [trial_left] = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        trial_right = trial_left;
        disp('Unsupported stimulus type for PlayFile test; reset to BROADBAND');
end

% remove any DCoffset
trial_left = trial_left - mom(trial_left,1);
trial_right = trial_right - mom(trial_right,1);

% normalize to ACPower (added 3/5/07 as in 2-source)
if ~XStimParams.HiDynamicRange & scaleFlag
    trial_left = trial_left / mom(trial_left,2);
    trial_right = trial_right / mom(trial_right,2);
end

% modulate stim1
if get(H.stim_type,'Value') ~= 9 & ~strcmp(XStimParams.mod_type,'None')
    Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
    trial_right = trial_right .* Envelope(:)';
    trial_left = trial_left .* Envelope(:)';
end

if get(H.stim_type,'Value') ~= 9
    %Ramp the stimuli
    ramp_time = 5; %ms
    [trial_left] = ramp_sound(trial_left,TDT.Fs,ramp_time);
    [trial_right] = ramp_sound(trial_right,TDT.Fs,ramp_time);
end

% remove any DCoffset
trial_left = trial_left - mom(trial_left,1);
trial_right = trial_right - mom(trial_right,1);

%Apply ITD filtering 
if (ITD ~=0)
    itdleft = 0; itdright = 0;
    ITD = round(str2num(get(H.PlayFile_ITD,'String')));
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
    trial_left = conv(trial_left,itd_filt_left);
    trial_right = conv(trial_right,itd_filt_right);
end

% Apply ILD
if(ILD ~=0)
    %%% DO THIS
end

% adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
% presented through DSPs (*.eq) at 0,0 (broadband)
trial_left = trial_left * TDT.hrtf_Lfactor;
trial_right = trial_right * TDT.hrtf_Rfactor;

%scale stimuli 3/7/07
% note scale factor changes with XStimParams.HiDynamicRange
if scaleFlag
trial_left = trial_left * scaleFactor;
trial_right = trial_right * scaleFactor;
end

%Add in the leading silent period
trial_left =  [zeros(1,silence_len) trial_left];
trial_right = [zeros(1,silence_len) trial_right];

%Add in the trailing silent period
trial_left =  [trial_left zeros(1,silence_len2)];
trial_right = [trial_right zeros(1,silence_len2)];

% save stimuli to disk
S232('push16',trial_left,length(trial_left));
S232('qpop16',BUF.L1);
fname = ['PlayFile1_' str1 '.left_1'];
evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
        [' ''' FN.temp_stim_path fname ''' '] ',0);'];
eval(evalstr);
%temp_left = dama2pc(BUF.L1);
S232('push16',trial_right,length(trial_right));
S232('qpop16',BUF.R1);
fname = ['PlayFile1_' str1 '.right_1'];
evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
        [' ''' FN.temp_stim_path fname ''' '] ',0);'];
eval(evalstr);
%temp_right = dama2pc(BUF.R1);

% set up for next trial                        
set(H.PlayFile_remtrials,'String','0');
set(H.PlayFile_status,'BackgroundColor','blue');
pause(0);

%%%%%%%%%%%%%%%%%%% finished making stimuli
%Begin playing sounds   
set(H.PlayFile_status,'String','Status: playing stimuli');
set(H.PlayFile_status,'BackgroundColor','green');
set(H.PlayFile_remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
buffcycle = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.PlayFilefig') & get(H.PlayFile_recorddata,'Value'))
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
%loop for reps
while (exist1('H.PlayFilefig') & (repnum <= XStimParams.numreps))
    spikes_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    tic
    %Check for pause by user
    if pause_check    return; end
    
    % load stimuli
    fname = ['PlayFile1_' str1 '.left_1'];
    evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
    eval(evalstr);
    evalstr = ['S232(''qpop16'''  ',BUF.L1);'];
    eval(evalstr);
    fname = ['PlayFile1_' str1 '.right_1'];
    evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
    eval(evalstr);
    evalstr = ['S232(''qpop16'''  ',BUF.R1);'];
    eval(evalstr);
    
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
    
    %Check for pause by user
    pause(DUR/1000+.2);     % to bypass PD1status bug
    while(S232('PD1status',1)) usec_delay(1000);  end
    S232('PD1stop',1);
    
    %Stop the m110 and get spikes
    m110dx( C_.STOP);
    spikes = m110dx( C_.DATA, XStimParams.curr_stimdur); 			% Take XStimParams.curr_stimdur spikes max
    ind = find(spikes ~= 0); 						% Get clock events that are spikes
    spikes = spikes(ind);
    ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
    if(exist1('H.PlayFilefig') & ~isempty(spikes)) 
        spikes = spikes(ind);
        spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
        repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
        Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
        diff(spikes_trial)
    end
    finalspikematrix = finalspikematrix + length(spikes);
    if pause_check    return; end
    
    set(H.PlayFile_remtrials,'String','1');
    pause(0);
    
    %Record Data
    if(exist1('H.PlayFilefig') & get(H.PlayFile_recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial]];
        record_data3(XStimParams,datamatrix);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.PlayFile_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure

set(H.PlayFile_status,'String',['Status: results: ' num2str(finalspikematrix) ' spikes']);
set(H.PlayFile_status,'BackgroundColor','blue');

set(H.exitPlayFile,'Visible','on');
set(H.resetPlayFile,'Visible','on');

% increment test number
if(exist1('H.PlayFilefig') & get(H.PlayFile_recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end


%%%%%%%%%%%
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
    case 'Sq wave'
        T = 0:Inc:(DUR/1000 - Inc);
        Tone = (param1 / 2)* square(2 * pi * param2 .* T + (.75 * 2 * pi + param3));
        Envelope = Tone + (1-param1/2);
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
                    set(H.PlayFile_modfile,'String',[FN.mod_path FN.mod]);
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
while (exist1('H.PlayFilefig') & get(H.pausePlayFile,'Value'))
    pause(0);
    if(~exist1('H.PlayFilefig')) return; end         
    set(H.exitPlayFile,'Visible','on');
    set(H.resetPlayFile,'Visible','on');
    if(exist1('H.PlayFilefig') & get(H.resetPlayFile,'Value') == 1)
        set(H.resetPlayFile,'Value',0);
        set(H.pausePlayFile,'Value',0);
        Reset_PlayFile;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_PlayFile;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end