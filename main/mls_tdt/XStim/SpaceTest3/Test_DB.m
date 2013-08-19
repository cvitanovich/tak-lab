function [] = test_DB()

globals_mii
globals_tdt
load e:\kip\matlab\mls_tdt\XStim\control\FN_current
load e:\kip\matlab\mls_tdt\XStim\control\XStimParams_current

% set params
XStimParams.curr_freq = 3000;
XStimParams.curr_ABL = 40;
XStimParams.curr_stimdur = 3000;
XStimParams.numreps = 4;
stim_type = 1;      % 1) tone  2) BBN
dec_factor = 2;


% get HRTFs
dir = 0;
eval(['load -mat ' FN.space_path FN.space_eq]);
locind = max(find(dir(1,:) == 0 & dir(2,:) == 0));
HRTF_L = TF1(locind,:)';
HRTF_R = TF2(locind,:)';
clear dir TF1 TF2

% Initialize Application and get AP2 and XBUS locks
if(S232('S2init', 0, 'INIT_PRIMARY', 20000) == 0)
   disp('Cannot initialize a secondary process')
   return;
end

if(S232('APlock', 100, 0) == 0)
  disp('Cannot acquire lock on AP2 Card')
  s232('S2close');
  return;
end

if(S232('XBlock', 100, 0) == 0)
  disp('Cannot acquire lock on X Bus')
  s232('APunlock', 0);
  s232('S2close');
  return;
end


% setup PD1
S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);

% reset DSPs
S232('PD1nstrms',1,2,2);		% 2 input and 2 output buffers
S232('PD1resetDSP',1,hex2dec('FFF'));
S232('dropall');
S232('PD1clrsched',1);

% for left ear: IB[0] >> IREG[0] >> DSP[0] >>  DAC[0]
S232('PD1addsimp',1,S232('DSPout',0),S232('DAC',0));
S232('PD1addsimp',1,S232('IREG',0),S232('DSPin',0));  
S232('PD1specIB',1,S232('IB',0),S232('IREG',0)); %IB to DSPin

% for right ear:IB[1] >> IREG[0] >> DSP[1] >>  DAC[1]
S232('PD1addsimp',1,S232('DSPout',1),S232('DAC',1));
S232('PD1addsimp',1,S232('IREG',1),S232('DSPin',1));  
S232('PD1specIB',1,S232('IB',1),S232('IREG',1));

% setup 2x OBs to ADC
s232('PD1specOB', 1,S232('OB',0),S232('ADC',0));
s232('PD1specOB', 1,S232('OB',1),S232('ADC',1));

%Specify DAMA buffers
clear BUF
BUF.L1			= 1;
BUF.R1			= 2;
BUF.L2			= 3;
BUF.R2			= 4;
BUF.playseq_L	= 5;
BUF.playseq_R	= 6;
BUF.playspec	= 7;
S232('allot16',BUF.L1, TDT.buf_pts);
S232('allot16',BUF.R1, TDT.buf_pts);
S232('allot16',BUF.L2, TDT.buf_pts);
S232('allot16',BUF.R2, TDT.buf_pts);

BUF.recL1		= 11;
BUF.recR1		= 12;
BUF.recL2		= 13;
BUF.recR2		= 14;
BUF.recseq_L	= 15;
BUF.recseq_R	= 16;
BUF.recspec		= 17;
S232('allot16',BUF.recL1, TDT.buf_pts);
S232('allot16',BUF.recR1, TDT.buf_pts);
S232('allot16',BUF.recL2, TDT.buf_pts);
S232('allot16',BUF.recR2, TDT.buf_pts);

BUF.dec1        = 21;
BUF.dec2        = 22;
S232('allot16',BUF.dec1, ceil(TDT.buf_pts/(2^dec_factor)));
S232('allot16',BUF.dec2, ceil(TDT.buf_pts/(2^dec_factor)));

%Make PLAY sequence buffers
S232('dropall');
S232('allot16',BUF.playseq_L,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.L1);
S232('make',1,1);
S232('make',2,BUF.L2);
S232('make',3,1);
S232('make',4,0);
S232('qpop16',BUF.playseq_L);

S232('allot16',BUF.playseq_R,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.R1);
S232('make',1,1);
S232('make',2,BUF.R2);
S232('make',3,1);
S232('make',4,0);
S232('qpop16',BUF.playseq_R);

%Make play SPECification buffer
S232('allot16',BUF.playspec,11);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L);
S232('make',1,BUF.playseq_R);
S232('make',2,0);
S232('qpop16',BUF.playspec);

%Make RECORD sequence buffers
S232('dropall');
S232('allot16',BUF.recseq_L,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.recL1);
S232('make',1,1);
S232('make',2,BUF.recL2);
S232('make',3,1);
S232('make',4,0);
S232('qpop16',BUF.recseq_L);

S232('allot16',BUF.recseq_R,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.recR1);
S232('make',1,1);
S232('make',2,BUF.recR2);
S232('make',3,1);
S232('make',4,0);
S232('qpop16',BUF.recseq_R);

%Make play SPECification buffer
S232('allot16',BUF.recspec,11);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.recseq_L);
S232('make',1,BUF.recseq_R);
S232('make',2,0);
S232('qpop16',BUF.recspec);

%Set MII parameters
mii_us_per_sample = 10; 				%microsecond per sample
mii_separation = 100; 					%only take events separated by 100 samples (i.e., 1 ms)
s_pt=0;

% set double buffering params
recordBUFIDs = [BUF.recL1 BUF.recL2 BUF.recR1 BUF.recR2];
decimationBUFID_factors = [BUF.dec1 BUF.dec2 dec_factor dec_factor];
playBUFIDs = [BUF.L1 BUF.L2 BUF.R1 BUF.R2];
bufSTARTPOSITIONs = [0 0 (s_pt + TDT.buf_pts) (s_pt + TDT.buf_pts)];
bufLENGTHs = [0 0 TDT.buf_pts TDT.buf_pts];             % if value==0, then don't use this BUF

% set ABL
%S232('PA4mute',1);
%S232('PA4mute',2);

S232('PA4atten',1,abs(XStimParams.curr_ABL)-20);
S232('PA4atten',2,abs(XStimParams.curr_ABL)-20);

silence_len = (XStimParams.silence_lead * round(TDT.Fs/1000));
silence_len2 = (XStimParams.silence_trail * round(TDT.Fs/1000));

if 0    %%%%%% Make stimulus
    switch stim_type    
        case 1      % tone
            sound_left = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur,0);
        case 8        %Broadband Noise
            sound_left = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        otherwise
            disp('wrong stimtype')
    end
    % remove DC
    sound_left = sound_left(:)' - mean(sound_left);
    % scale to ACPower = 1;
    sound_left = sound_left / mom(sound_left,2);
    %reset stim-dur for actual sound made
    XStimParams.curr_stimdur = length(sound_left)* (1000/TDT.Fs);
    % Ramp the stimulus
    sound_left = ramp_sound(sound_left,TDT.Fs,5);
    %Add silent period prior to stimulus to calculate spontaneous rate
    sound_left =  [zeros(1,silence_len) sound_left];
    %Add silent period after stimulus for HRTF filtering
    sound_left =  [sound_left zeros(1,TDT.hrtffiltlen)];
    %Add silent period after stimulus for additional spike gathering
    sound_left =  [sound_left zeros(1,silence_len2)];
    %scale stimulus
    sound_left = sound_left * TDT.scaleFactor;
    % R=L
    sound_right = sound_left;
end

% calc npts to play
npts_totalplay = silence_len + (XStimParams.curr_stimdur * round(TDT.Fs/1000))+ TDT.hrtffiltlen + silence_len2;

S232('PD1npts',1,npts_totalplay);	
S232('dropall');
% preload HRTFs to DSPs nearest to DACs
S232('PD1resetDSP',1, hex2dec('1'));
S232('pushf',HRTF_L(:,1)',255);
S232('PreLoadRaw',1,S232('DSPid',0),'MONO','STACK','','',TDT.hrtf_scale,TDT.hrtf_scale,1);
S232('PD1resetDSP',1, hex2dec('2'));
S232('pushf',HRTF_R(:,1)',255);
S232('PreLoadRaw',1,S232('DSPid',1),'MONO','STACK','','',TDT.hrtf_scale,TDT.hrtf_scale,1);


%%%%%%%%%%%%%%%%%%% play one blank sound first to work-around TDT glitch
S232('seqplay',BUF.playspec);
%Set up MII
m100x( C_.INIT );
m110dx( C_.INIT );
m110dx( C_.CLOCK, mii_us_per_sample);
m110dx( C_.MODE, M_.PST );

%Fill first stimulus buffers
if 0
s_pt = 1;
S232('disk2dama16', BUF.L1,FNtempL,s_pt);
S232('disk2dama16', BUF.R1,FNtempR,s_pt);
%Fill second stimulus buffers
s_pt = (s_pt + TDT.buf_pts);
S232('disk2dama16', BUF.L1,FNtempL,s_pt);
S232('disk2dama16', BUF.R1,FNtempR,s_pt);
end

% fill buffers
S232('dropall');
S232('dpush',TDT.buf_pts);
S232('tone',XStimParams.curr_freq,1e6/TDT.Fs);
S232('scale',10000);
S232('qdup');
S232('qpop16',BUF.L1);
S232('qdup');
S232('qpop16',BUF.R1);

S232('tone',XStimParams.curr_freq/2,1e6/TDT.Fs);
S232('scale',10000);
S232('qdup');
S232('qpop16',BUF.L2);
S232('qpop16',BUF.R2);
S232('dropall');

%%% arming, firing and triggering from MII port0 all from within
%%% double_buffer.dll
double_buffer('e:\kip\noisetc','e:\kip\noisetc', 'e:\kip\noisetc','e:\kip\noisetc', ...
    recordBUFIDs,decimationBUFID_factors,playBUFIDs, ...
    bufSTARTPOSITIONs,bufLENGTHs,npts_totalplay);    

%Stop the m110 and get spikes
m110dx( C_.STOP);
spikes = m110dx( C_.DATA, 1000); 			%Take 1000 spikes max

pause(1)

%%%%% set up for looping
S232('PA4atten',1,abs(XStimParams.curr_ABL)-20);
S232('PA4atten',2,abs(XStimParams.curr_ABL)-20);

repnum = 1;

%%% start looping through reps 
while (repnum <= XStimParams.numreps)
    
    %%% start looping through trials (locations)
    S232('dropall');
    
    % fill buffers
    S232('dropall');
    S232('dpush',TDT.buf_pts);
    S232('tone',XStimParams.curr_freq,1e6/TDT.Fs);
    S232('scale',10000);
    S232('qdup');
    S232('qpop16',BUF.L1);
    S232('qdup');
    S232('qpop16',BUF.R1);
    
    S232('tone',XStimParams.curr_freq/2,1e6/TDT.Fs);
    S232('scale',10000);
    S232('qdup');
    S232('qpop16',BUF.L2);
    S232('qpop16',BUF.R2);
    S232('dropall');
    
    double_buffer('e:\kip\noisetc','e:\kip\noisetc', 'e:\kip\noisetc','e:\kip\noisetc', ...
        recordBUFIDs,decimationBUFID_factors,playBUFIDs, ...
        bufSTARTPOSITIONs,bufLENGTHs,npts_totalplay);    
    
    %Stop the m110 and get spikes
    m110dx( C_.STOP);
    repnum = repnum + 1;
    
    
    pause(5)
end


if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('PA4mute',1);
   S232('PA4mute',2);
   S232('PD1stop',1);
   S232('PD1clear',1);
end
m110dx( C_.STOP);

if(S232('getAPlockstatus') & S232('getXBlockstatus'))
   S232('trash');
   S232('dropall');
   S232('S2close');
end
