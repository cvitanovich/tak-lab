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
    [HRTF.left,HRTF.right] = readHRTFdotMAT(PDR.HRTF_fname,el,az);
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
tmp = makeTest(PDR.TEST_seed,1000,bandwidth(1),bandwidth(2),PDR.stim_Fs,PDR.TEST_ramp,PDR.TEST_base_rms);
tmp = (0.999)*tmp ./ (max(abs(tmp))); % NORMALIZE
% Convolve with HRTF (ABL Equalized) for frontal location (El,Az)=(0,0)
PDR.OCTAVE_left = filter(HRTF.left,1,tmp);
PDR.OCTAVE_right = filter(HRTF.right,1,tmp);

%% MAKE BROADBAND NOISE (2-11kHz)
bandwidth=[2000 11000];
tmp= makeTest(PDR.TEST_seed,1000,bandwidth(1),bandwidth(2),PDR.stim_Fs,PDR.TEST_ramp,PDR.TEST_base_rms);
tmp = (0.999)*tmp ./ (max(abs(tmp))); % NORMALIZE
% Convolve with HRTF (ABL Equalized) for frontal location (El,Az)=(0,0)
PDR.BBN_left = filter(HRTF.left,1,tmp);
PDR.BBN_right = filter(HRTF.right,1,tmp);

%% MAKE GAMMATONE (cF=6kHz)
PDR.GTONE_cF = 6; PDR.GTONE_species = 'owl';
PDR.GTONE_coefs = makeGammaFIR(PDR.stim_Fs,PDR.GTONE_cF,PDR.GTONE_species);
tmp=rand(1,PDR.buf_pts);
tmp=filter(PDR.GTONE_coefs,1,tmp);
% Convolve with HRTF (ABL Equalized) for frontal location (El,Az)=(0,0)
PDR.GTONE_left = filter(HRTF.left,1,tmp);
PDR.GTONE_right = filter(HRTF.right,1,tmp);