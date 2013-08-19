function HRTF = make_calib_sounds
%% MAKE SOUND TOKENS FOR SPATIAL PDR CALIBRATION
global PDR

HRTF.directory=PDR.HRTF_directory;
HRTF.fname=PDR.HRTF_fname;
% preallocate a sound buffer:
tmp=NaN*ones(1,PDR.buf_pts);

%% LOAD HRTF COEFFICIENTS
if(strcmp(PDR.HRTF_fname((end-3):end),'.mat'))
    el=0; az=0;
    [HRTF.left,HRTF.right] = readHRTFdotMAT(PDR.HRTF_directory,PDR.HRTF_fname,el,az);
else % using other format (e.g. for 930 or 929)
    %read HRTF coefficients files:
    
    HRTF = MTLreadHDR(HRTF);
    HRTF = MTLreadDIR(HRTF);
    % THIS CONVERSION IS NOT NECESSARY WITH ANY OF THE BIRDS WE USE NOW:
    % %convert coordinates to double polar, only if not using 929 or 930 or ones
    % if (~strcmp(PDR.HRTF_fname(1:3),'930') & ~strcmp(PDR.HRTF_fname(1:3),'929'))
    %     PDR.HRTF_dir_matrix = sphere2double(PDR.HRTF_dir_matrix);
    % end
    
    HRTF.left = nan*ones(HRTF.nlines,1);
    HRTF.right = HRTF.left;
    direc = HRTF.dir_matrix; % NOTE: first row = Elevation and 2nd row = Azimuth !!!
    [r c]=find(direc(1,:)==0);
    idx=find(direc(2,c)==0); % using frontal location for calibration
    HRTF.left = MTLreadCH(idx*2-1, HRTF);
    HRTF.right = MTLreadCH(idx*2, HRTF);
end

%% MAKE NARROWBAND (Octave Band: 4-8kHz) NOISE
bandwidth=[4000 8000];
ok=0; cnt=0;
while(~ok)
    PDR.OCTAVE_seed=PDR.TEST_seed+cnt;
    tmp= makeTest(PDR.OCTAVE_seed,PDR.TEST_dur,bandwidth(1),bandwidth(2),PDR.stim_Fs,0);
    tmp = 0.999.*(tmp ./ (max(abs(tmp)))); % NORMALIZE
    % Convolve with HRTF (ABL Equalized) for frontal location (El,Az)=(0,0)
    PDR.OCTAVE_left = filter(HRTF.left,1,tmp);
    PDR.OCTAVE_right = filter(HRTF.right,1,tmp);
    % NORMALIZE TO +/- 1
    mx=max(abs([PDR.OCTAVE_left PDR.OCTAVE_right]));
    PDR.OCTAVE_left=PDR.OCTAVE_left./mx;
    PDR.OCTAVE_right=PDR.OCTAVE_right./mx;
    rms_left = sqrt(mean(PDR.OCTAVE_left.^2));
    rms_right = sqrt(mean(PDR.OCTAVE_right.^2));
    rms_avg = (rms_left+rms_right)/2;
    
    PDR.OCTAVE_left=(PDR.TEST_target_rms/rms_avg).*PDR.OCTAVE_left; % sets rms ABL to target rms
    PDR.OCTAVE_right=(PDR.TEST_target_rms/rms_avg).*PDR.OCTAVE_right; % sets rms ABL to target rms
    if(max(abs([PDR.OCTAVE_left PDR.OCTAVE_right]))<1)
        ok=1; % avoids CLIPPING!
        break
    end
    cnt=cnt+1;
end
PDR.OCTAVE_left = rampMySound(PDR.OCTAVE_left,PDR.TEST_ramp,PDR.stim_Fs);
PDR.OCTAVE_right = rampMySound(PDR.OCTAVE_right,PDR.TEST_ramp,PDR.stim_Fs);

%% MAKE BROADBAND NOISE (2-11kHz)
bandwidth=[2000 11000];
ok=0; cnt=0;
while(~ok)
    PDR.BBN_seed=PDR.TEST_seed+cnt;
    tmp= makeTest(PDR.BBN_seed,PDR.TEST_dur,bandwidth(1),bandwidth(2),PDR.stim_Fs,0);
    tmp = 0.999.*(tmp ./ (max(abs(tmp)))); % NORMALIZE
    % Convolve with HRTF (ABL Equalized) for frontal location (El,Az)=(0,0)
    PDR.BBN_left = filter(HRTF.left,1,tmp);
    PDR.BBN_right = filter(HRTF.right,1,tmp);
    % NORMALIZE TO +/- 1
    mx=max(abs([PDR.BBN_left PDR.BBN_right]));
    PDR.BBN_left= PDR.BBN_left./mx;
    PDR.BBN_right= PDR.BBN_right./mx;
    rms_left= sqrt(mean(PDR.BBN_left.^2));
    rms_right= sqrt(mean(PDR.BBN_right.^2));
    rms_avg = (rms_left+rms_right)/2;
    PDR.BBN_left=(PDR.TEST_target_rms/rms_avg).*PDR.BBN_left; % sets rms ABL to target rms
    PDR.BBN_right=(PDR.TEST_target_rms/rms_avg).*PDR.BBN_right; % sets rms ABL to target rms
    if(max(abs([PDR.BBN_left PDR.BBN_right]))<1)
        ok=1; % avoids CLIPPING!
        break
    end
    cnt=cnt+1;
end

PDR.BBN_left = rampMySound(PDR.BBN_left,PDR.TEST_ramp,PDR.stim_Fs);
PDR.BBN_right = rampMySound(PDR.BBN_right,PDR.TEST_ramp,PDR.stim_Fs);

%% MAKE GAMMATONE (cF=6kHz)
PDR.GTONE_cF = 6000; PDR.GTONE_species = 'owl';
PDR.GTONE_coefs = makeGammaFIR(PDR.stim_Fs,PDR.GTONE_cF,PDR.GTONE_species);
ok=0; cnt=0;
while(~ok)
    PDR.GTONE_seed=PDR.TEST_seed+cnt;
    rand('state',PDR.GTONE_seed);
    tmp=rand(1,PDR.buf_pts);
    tmp=filtfilt(PDR.GTONE_coefs,1,tmp);
    tmp = 0.999.*(tmp ./ (max(abs(tmp)))); % NORMALIZE
    % Convolve with HRTFs (ABL Equalized) for frontal location (El,Az)=(0,0)
    PDR.GTONE_left = filter(HRTF.left,1,tmp);
    PDR.GTONE_right = filter(HRTF.right,1,tmp);
    % NORMALIZE TO +/- 1
    mx=max(abs([PDR.GTONE_left PDR.GTONE_right]));
    PDR.GTONE_left= PDR.GTONE_left./mx;
    PDR.GTONE_right= PDR.GTONE_right./mx;
    rms_left= sqrt(mean(PDR.GTONE_left.^2));
    rms_right= sqrt(mean(PDR.GTONE_right.^2));
    rms_avg = (rms_left+rms_right)/2;
    PDR.GTONE_left=(PDR.TEST_target_rms/rms_avg).*PDR.GTONE_left; % sets rms ABL to target rms
    PDR.GTONE_right=(PDR.TEST_target_rms/rms_avg).*PDR.GTONE_right; % sets rms ABL to target rms
    if(max(abs([PDR.GTONE_left PDR.GTONE_right]))<1)
        ok=1; % avoids CLIPPING!
        break
    end
    cnt=cnt+1;
end
PDR.GTONE_left = rampMySound(PDR.GTONE_left,PDR.TEST_ramp,PDR.stim_Fs);
PDR.GTONE_right = rampMySound(PDR.GTONE_right,PDR.TEST_ramp,PDR.stim_Fs);

