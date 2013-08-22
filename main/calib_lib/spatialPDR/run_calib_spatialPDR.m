%% main script to calibrate for spatial PDR (depends on calib_sound.m)
%% should calibrate knowles mics and spl meter beforehand

%% PARAMETERS:
% note: microphones should also be calibrated for SPL vs. Voltage on the same day
C.data_path='C:\andrew\data\calib\';
% attens to try
C.atten_list=25:-5:5;
% what scales to try?
C.scales_2_try_for_cutoffs=fliplr(10.^([0:.09:4.5 log10(32760)]));
C.nscales=50;
C.stim_Fs=30000;
C.buf_pts=30000; % 1 sec sounds
C.stim_dur = C.buf_pts/C.stim_Fs;
% HRTF coefficients file:
C.HRTF_fname='1073AC_eq_ABLequal.mat'; % what file to use for HRTFs???
C.HRTF_nlines=255;
C.HRTF_directory='C:\Documents and Settings\andrew\My Documents\GitHub\tak-lab\main\HRTFs\Matlab_V6\';
C.HRTF_Left=[]; % left coefs
C.HRTF_Right=[]; % right coefs

% sounds
C.sounds = {'GTONE','OCTAVE','BBN'};
C.ramp = 5; % 5 ms ramp
C.target_rms=0.1; % desired rms amplitude of sound vector
C.GTONE_cF=6000;
C.GTONE_species = 'owl';
C.OCTAVE_bandwidth=[4000 8000];
C.BBN_bandwidth=[2000 11000];
% filename and data path
c=clock; c=round(c);
yr=num2str(c(1)); yr=yr(end-1:end);
mo=num2str(c(2)); day=num2str(c(3));


%% LOAD KNOWLES VOLTAGE CALIBRATIONS:
uiopen([C.data_path 'knowles*.mat']);
C.knowles_file = [KNOWLES.data_path KNOWLES.filename];
disp(['Loaded Knowles Mic calibrations from: ' C.knowles_file]);

%% LOAD SPL METER VOLTAGE CALIBRATIONS:
uiopen([C.data_path 'spl_meter*.mat']);
C.spl_meter_file = [SPL_METER.data_path SPL_METER.filename];
disp(['Loaded SPL Meter calibrations from: ' C.spl_meter_file]);

%% LOAD HRTF COEFFICIENTS:
if(strcmp(C.HRTF_fname((end-3):end),'.mat')) % DOT MAT FORMAT
    LT=zeros(1,C.HRTF_nlines);
    RT=LT;
    EL=0; AZ=0;
    [LT, RT] = readHRTFdotMAT(C.HRTF_directory,C.HRTF_fname,EL,AZ);
    HRTF.left=LT;
    HRTF.right=RT;
else % using another format (e.g. for 930 or 929)
    %read HRTF coefficients files:
    [HRTF] = MTLreadHDR(C.HRTF_directory,C.HRTF_fname);
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
    
    % HRTF filtered Version:
    % Convolve with HRTF (ABL Equalized) for frontal location (El,Az)=(0,0)
    C.OCTAVE_left = filter(HRTF.left,1,tmp);
    C.OCTAVE_right = filter(HRTF.right,1,tmp);
    % NORMALIZE TO +/- 1
    mx=max(abs([C.OCTAVE_left C.OCTAVE_right]));
    C.OCTAVE_left=C.OCTAVE_left./mx;
    C.OCTAVE_right=C.OCTAVE_right./mx;
    rms_left = sqrt(mean(C.OCTAVE_left.^2));
    rms_right = sqrt(mean(C.OCTAVE_right.^2));
    rms_avg = (rms_left+rms_right)/2;
    C.OCTAVE_left=(C.target_rms/rms_avg).*C.OCTAVE_left; % sets rms ABL to target rms
    C.OCTAVE_right=(C.target_rms/rms_avg).*C.OCTAVE_right; % sets rms ABL to target rms
    
    % Version w/o HRTF Filtering:
    C.OCTAVE_sound = tmp;
    rms = sqrt(mean(C.OCTAVE_sound.^2));
    C.OCTAVE_sound = (C.target_rms/rms).*C.OCTAVE_sound; % set rms equal to target rms
    if(max(abs(C.OCTAVE_sound))<1 && max(abs([C.OCTAVE_left C.OCTAVE_right]))<1)
        ok=1; % avoids CLIPPING!
        break
    end
    cnt=cnt+1;
end
C.OCTAVE_left = rampMySound(C.OCTAVE_left,C.ramp,C.stim_Fs);
C.OCTAVE_right = rampMySound(C.OCTAVE_right,C.ramp,C.stim_Fs);
C.OCTAVE_sound = rampMySound(C.OCTAVE_sound,C.ramp,C.stim_Fs);

%% MAKE BROADBAND NOISE
bandwidth=C.BBN_bandwidth;
ok=0; cnt=0; seed=1;
while(~ok)
    C.BBN_seed=seed+cnt;
    tmp= makeTest(C.BBN_seed,C.stim_dur*1000,bandwidth(1),bandwidth(2),C.stim_Fs,0);
    tmp = 0.999.*(tmp ./ (max(abs(tmp)))); % NORMALIZE
    
    % HRTF Version:
        % Convolve with HRTF (ABL Equalized) for frontal location (El,Az)=(0,0)
        C.BBN_left = filter(HRTF.left,1,tmp);
        C.BBN_right = filter(HRTF.right,1,tmp);
        % NORMALIZE TO +/- 1
        mx=max(abs([C.BBN_left C.BBN_right]));
        C.BBN_left= C.BBN_left./mx;
        C.BBN_right= C.BBN_right./mx;
        rms_left= sqrt(mean(C.BBN_left.^2));
        rms_right= sqrt(mean(C.BBN_right.^2));
        rms_avg = (rms_left+rms_right)/2;
        C.BBN_left=(C.target_rms/rms_avg).*C.BBN_left; % sets rms ABL to target rms
        C.BBN_right=(C.target_rms/rms_avg).*C.BBN_right; % sets rms ABL to target rms
    
    % Version without HRTF Filtering
    C.BBN_sound = tmp;
    rms = sqrt(mean(C.BBN_sound.^2));
    C.BBN_sound = (C.target_rms/rms).*C.BBN_sound; % set rms equal to target
    
    if(max(abs(C.BBN_sound)<1) && max(abs([C.BBN_left C.BBN_right]))<1)
        ok=1; % avoids CLIPPING!
        break
    end
    
    cnt=cnt+1;
end
C.BBN_left = rampMySound(C.BBN_left,C.ramp,C.stim_Fs);
C.BBN_right = rampMySound(C.BBN_right,C.ramp,C.stim_Fs);
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
    
    % HRTF Filtered Version
    % Convolve with HRTFs (ABL Equalized) for frontal location (El,Az)=(0,0)
    C.GTONE_left = filter(HRTF.left,1,tmp);
    C.GTONE_right = filter(HRTF.right,1,tmp);
    % NORMALIZE TO +/- 1
    mx=max(abs([C.GTONE_left C.GTONE_right]));
    C.GTONE_left= C.GTONE_left./mx;
    C.GTONE_right= C.GTONE_right./mx;
    rms_left= sqrt(mean(C.GTONE_left.^2));
    rms_right= sqrt(mean(C.GTONE_right.^2));
    rms_avg = (rms_left+rms_right)/2;
    C.GTONE_left=(C.target_rms/rms_avg).*C.GTONE_left; % sets rms ABL to target rms
    C.GTONE_right=(C.target_rms/rms_avg).*C.GTONE_right; % sets rms ABL to target rms

    % Version w/o HRTF Filtering
    C.GTONE_sound=tmp;
    rms=sqrt(mean(C.GTONE_sound.^2));
    C.GTONE_sound = (C.target_rms/rms).*C.GTONE_sound; % set rms equal to target
    
    if(max(abs(C.BBN_sound)<1) && max(abs([C.GTONE_left C.GTONE_right]))<1)
        ok=1; % avoids CLIPPING!
        break
    end
    
    cnt=cnt+1;
end
C.GTONE_left = rampMySound(C.GTONE_left,C.ramp,C.stim_Fs);
C.GTONE_right = rampMySound(C.GTONE_right,C.ramp,C.stim_Fs);
C.GTONE_sound = rampMySound(C.GTONE_sound,C.ramp, C.stim_Fs);


%% NO BIRD, JUST SPEAKER AND SPL METER
C.HRTFs=0; C.EPHONES=0; C.SPEAKERS=1;
% mic_type = Knowles or SPL meter
C.MIC_TYPE='SPL Meter';
C.CALIB_TYPE='Intraural';
MIC.coeffs = SPL_METER.coeffs;
for j=length(C.sounds)
    % sound type
    snd=C.sounds{j};
    for k=length(C.atten_list)
        % atten
        C.current_atten=C.atten_list(k);
        % filenames
        C.filename=['Calib' C.CALIB_TYPE C.MIC_TYPE '_' snd '_' yr mo day '_Atten' num2str(C.current_atten)];
        % run calibration code
        C=calib_sound(C,MIC);
        % done. now mov on to the next atten/sound
    end
end

%% CALIBRATE INTRAURALLY, USING SPEAKERS
C.HRTFs=0; C.EPHONES=0; C.SPEAKERS=1;
% mic_type = Knowles or SPL meter
C.MIC_TYPE='Knowles';
C.CALIB_TYPE='Intraural';
MIC.coeffs=KNOWLES.coeffs;
for j=length(C.sounds)
    % sound type
    snd=C.sounds{j};
    for k=length(C.atten_list)
        % atten
        C.current_atten=C.atten_list(k);
        % filenames
        C.filename=['Calib' C.CALIB_TYPE C.MIC_TYPE '_' snd '_' yr mo day '_Atten' num2str(C.current_atten)];
        % run calibration code
        C=calib_sound(C,MIC);
        % done. now mov on to the next atten/sound
    end
end


%% CALIBRATE INTRAURALLY, USING HEADPHONES & HRTFs
C.HRTFs=1; C.EPHONES=1; C.SPEAKERS=0;
% mic_type = Knowles or SPL meter
C.MIC_TYPE='Knowles';
C.CALIB_TYPE='Intraural';
MIC.coeffs=KNOWLES.coeffs;
for j=length(C.sounds)
    % sound type
    snd=C.sounds{j};
    for k=length(C.atten_list)
        % atten
        C.current_atten=C.atten_list(k);
        % filenames
        C.filename=['Calib' C.CALIB_TYPE C.MIC_TYPE '_' snd '_' yr mo day '_Atten' num2str(C.current_atten)];
        % run calibration code
        C=calib_sound(C,MIC);
        % done. now mov on to the next atten/sound
    end
end