%% main script to calibrate for spatial PDR (depends on calib_sound.m)
%% should calibrate knowles mics and spl meter beforehand

% refresh
clear all; close all;

%% PARAMETERS:
% note: microphones should also be calibrated for SPL vs. Voltage on the same day
C.data_path='E:\andrew\data\calib\';
C.scales=300:300:32700; % range of scales to try
C.stim_Fs=30000;
C.buf_pts=15000; % 0.5 sec sounds
C.stim_dur = C.buf_pts/C.stim_Fs;
C.bird_id = 925;

% HRTF coefficients file:
HRTF.fname='925AD_eq_ABLequal.mat'; % what file to use for HRTFs???
HRTF.nlines=255;
HRTF.directory='C:\andrew\code\tak-lab\main\HRTFs\Matlab_V6\';

% check bird/HRTF file
btn=questdlg(['Using ' HRTF.fname ' with bird #' num2str(C.bird_id) ', okay?'],'Checkpoint','YES','NO','YES');
if(strcmp(btn,'NO')); return; end;

% sounds
C.sounds = {'BBN','GTONE','OCTAVE'};
C.ramp = 5; % 5 ms ramp
C.target_rms=0.1; % desired rms amplitude of sound vector
C.GTONE_cF=6000;
C.GTONE_species = 'owl';
C.OCTAVE_bandwidth=[4000 8000];
C.BBN_bandwidth=[2000 11000];

% date
c=clock; c=round(c);
yr=num2str(c(1)); yr=yr(end-1:end);
mo=num2str(c(2)); day=num2str(c(3));

% make a folder for these calibrations
parent_dir=C.data_path;
new_dir = ['calib' yr mo day];
if(~exist([parent_dir new_dir],'dir'))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(parent_dir,new_dir);
else
    SUCCESS=1;
end
if(~SUCCESS)
    hWarn=warndlg(['FAILURE! ' MESSAGE]);
    return;
end
C.data_path=[parent_dir new_dir '\'];
disp(['Data will be saved in: ' C.data_path]);

%% LOAD HRTF COEFFICIENTS:
if(strcmp(HRTF.fname((end-3):end),'.mat')) % DOT MAT FORMAT
    LT=zeros(1,HRTF.nlines);
    RT=LT;
    EL=0; AZ=0;
    [LT, RT] = readHRTFdotMAT(HRTF.directory,HRTF.fname,EL,AZ);
    HRTF.left=LT;
    HRTF.right=RT;
else % using another format (e.g. for 930 or 929)
    %read HRTF coefficients files:
    [HRTF] = MTLreadHDR(HRTF.directory,HRTF.fname);
    [HRTF] = MTLreadDIR(HRTF);
    
    % THIS CONVERSION IS NOT NECESSARY WITH ANY OF THE BIRDS WE USE NOW:
    % %convert coordinates to double polar, only if not using 929 or 930 or ones
    % if (~strcmp(C.HRTF_fname(1:3),'930') & ~strcmp(C.HRTF_fname(1:3),'929'))
    %     C.HRTF_dir_matrix = sphere2double(C.HRTF_dir_matrix);
    % end
    
    direc = HRTF.dir_matrix; % NOTE: first row = Elevation and 2nd row = Azimuth !!!
    idx=find(direc(1,:)==0 & direc(2,:)==0);
    chan_left = MTLreadCH(idx*2-1, HRTF);
    chan_right = MTLreadCH(idx*2, HRTF);
    HRTF.left=chan_left;
    HRTF.right=chan_right;
end

%% MAKE OCTAVE BAND NOISE
bandwidth=C.OCTAVE_bandwidth;
ok=0; cnt=0; seed=1;
while(~ok)
    C.OCTAVE_seed=seed+cnt;
    tmp= makeTest(C.OCTAVE_seed,C.stim_dur*1000,bandwidth(1),bandwidth(2),C.stim_Fs,0);
    tmp = 0.999.*(tmp ./ (max(abs(tmp)))); % NORMALIZE
    C.OCTAVE_sound = tmp;
    rms = sqrt(mean(C.OCTAVE_sound.^2));
    C.OCTAVE_sound = (C.target_rms/rms).*C.OCTAVE_sound; % set rms equal to target rms
    if(max(abs(C.OCTAVE_sound))<1)
        ok=1; % avoids CLIPPING!
        break
    end
    cnt=cnt+1;
end
C.OCTAVE_sound = rampMySound(C.OCTAVE_sound,C.ramp,C.stim_Fs);

%% MAKE BROADBAND NOISE
bandwidth=C.BBN_bandwidth;
ok=0; cnt=0; seed=1;
while(~ok)
    C.BBN_seed=seed+cnt;
    tmp= makeTest(C.BBN_seed,C.stim_dur*1000,bandwidth(1),bandwidth(2),C.stim_Fs,0);
    tmp = 0.999.*(tmp ./ (max(abs(tmp)))); % NORMALIZE
    C.BBN_sound = tmp;
    rms = sqrt(mean(C.BBN_sound.^2));
    C.BBN_sound = (C.target_rms/rms).*C.BBN_sound; % set rms equal to target
    if(max(abs(C.BBN_sound)<1))
        ok=1; % avoids CLIPPING!
        break
    end
    cnt=cnt+1;
end
C.BBN_sound = rampMySound(C.BBN_sound,C.ramp,C.stim_Fs);

%% MAKE GAMMATONE
C.GTONE_coefs = makeGammaFIR(C.stim_Fs,C.GTONE_cF,C.GTONE_species);
ok=0; cnt=0; seed=1;
while(~ok)
    C.GTONE_seed=seed+cnt;
    rand('state',C.GTONE_seed);
    tmp=rand(1,C.buf_pts);
    tmp=filtfilt(C.GTONE_coefs,1,tmp);
    tmp = 0.999.*(tmp ./ (max(abs(tmp)))); % NORMALIZE
    C.GTONE_sound=tmp;
    rms=sqrt(mean(C.GTONE_sound.^2));
    C.GTONE_sound = (C.target_rms/rms).*C.GTONE_sound; % set rms equal to target
    if(max(abs(C.BBN_sound)<1))
        ok=1; % avoids CLIPPING!
        break
    end
    cnt=cnt+1;
end
C.GTONE_sound = rampMySound(C.GTONE_sound,C.ramp, C.stim_Fs);

CHOICE=questdlg('What is your wish?', ...
    'Select Calib Type', ...
    'Freefield/No Bird','Intraural/Speakers','Intraural/Earphones','Freefield/No Bird');

switch CHOICE
    case 'Freefield/No Bird'
        %% NO BIRD, JUST SPEAKER AND KNOWLES
        C.attens=10:5:40; % range of attenuations to try for speaker
        h=warndlg('Prepare for freefield calib with speaker.','warning'); uiwait(h);
        C.HRTFs=0; C.EPHONES=0; C.SPEAKERS=1; 
        C.BIRD=0; % flag for bird present
        % mic_type = Knowles or SPL meter
        C.MIC_TYPE='Knowles';
        C.CALIB_TYPE='Freefield';
        C=calib_sound(C);
    case 'Intraural/Speakers'
        %% CALIBRATE INTRAURALLY, USING SPEAKERS
        C.attens=10:5:40; % range of attenuations to try for speaker
        h=warndlg('Prepare for intraural calib with speaker.','warning'); uiwait(h);
        C.HRTFs=0; C.EPHONES=0; C.SPEAKERS=1; 
        C.BIRD=1; % flag for bird present
        % mic_type = Knowles or SPL meter
        C.MIC_TYPE='Knowles';
        C.CALIB_TYPE='Intraural';
        C=calib_sound(C);
    case 'Intraural/Earphones'
        %% CALIBRATE INTRAURALLY, USING HEADPHONES & HRTFs
        C.attens=20:5:50; % range of attenuations to try for ephones
        h=warndlg('Prepare for intraural calib with ephones.','warning'); uiwait(h);
        C.HRTFs=1; C.EPHONES=1; C.SPEAKERS=0;
        C.BIRD=1; % flag for bird present
        % mic_type = Knowles or SPL meter
        C.MIC_TYPE='Knowles';
        C.CALIB_TYPE='Intraural';
        C=calib_sound(C,HRTF);
    otherwise
        hWarn=warndlg('Invalid Option!');
end