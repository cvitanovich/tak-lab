% SETUP DEFAULTS FOR spatialPDR EXPERIMENT

% select folder location based on computer being used:
if(~ispc && strcmp(getenv('USER'),'cvitanovich'))
    code_path='/Users/cvitanovich/Documents/MATLAB/tak-lab/projects/spatial_adaptation/';
    data_path='/Users/cvitanovich/Documents/MATLAB/data/';
    hrtf_path='/Users/cvitanovich/Documents/MATLAB/tak-lab/HRTFs/Matlab_V6/';
    calib_path='/Users/cvitanovich/Documents/MATLAB/tak-lab/calib_data/spatial_PDR/';
else
    code_path='C:\andrew\CORE\tak-lab\projects\spatial_adaptation\';
    data_path='C:\andrew\pdr\data\';
    hrtf_path='C:\andrew\CORE\tak-lab\HRTFs\Matlab_V6\';
    calib_path='C:\andrew\CORE\tak-lab\calib_data\spatial_PDR\';
end
    
PDR = struct(...                            % MAIN PARAMETERS:
    'DEBUG',1,...
    'calib',0,...                           % flag for calibrations
    'virtual',1, ...                        % flag for virtual sound presentation
    'record',1, ...                        % flag for recording pupillometer output
    'flag_adapt',1,...                     % flag for adapting stimulus
    'ntrials',1400,...                       % # of trials in session
    'npretrials',5,...                     % # of trials before first test trial
    'n_test_trials',[],...                  % # of test trials 
    'buf_pts',16384,...                     % number of pts in each buffer (~0.5s)
    'buf_dur',[],...                        % buffer duration in ms
    'isi_buf',7,...                         % no. buffers between trials
    'isi_time',[],...                       % ISI (seconds) between trials
    'decimationfactor',4, ...               % decimation factor for data collection
    'stim_Fs', 30000, ...                   % sampling rate (Hz) for stimuli
    'npts_totalplay',[],...                % total no. pts in the session
    'len_session',[], ...                   % length of session (in minutes)
    'starttime',[], ...                     % session start time
    'stoptime', [], ...                     % session stop time
    'code_path', code_path,...  % path to code
    'data_path', data_path, ...     % data path (normally: c:\andrew\pdr\data normally, calib: c:\andrew\pdr\calib_data)
    'base_atten',60,...                      % base attenuation for PA4
    'filename',[],...                       % filename for session data
    'bird_id',925,...
    'hab_id',0,...
    ...
    ...                                     % ADAPTING STIMULUS PARAMETERS:
    'ADAPT_loc',[0 0],...                   % location of adapting stimulus
    'ADAPT_species','owl',...               % species
    'ADAPT_soundtype','GTONE',...            % type of sound for adapting stimulus
    'ADAPT_cF',6000,...                     % center frequency (Hz) (if relevant)
    'ADAPT_coefs',[],...                    % FIR coefficients to generate adapting stimulus
    'ADAPT_SPL',6,...                      % SPL for adaptor
    'ADAPT_scale',[],...                    % scale for adaptor (calculated using calibration data)
    'ADAPT_ramp',5,...                      % ramp for trial segments
    'ADAPT_nstates',16,...                  % no. of desired states for adaptor
    'ADAPT_dur',[],...                      % adaptor duration in seconds (calculated from isi_buf) in seconds
    'ADAPT_target_rms',0.1,...              % desired rms amplitude for adaptor (before scaling)
    'ADAPT_seeds',[],...                    % seed for generating adaptor
    ...
    ...                                     % TEST STIMULUS PARAMETERS:
    'TEST_soundtype','OCTAVE',...           % type of test sound ('gammatone', 'octave', 'broadband', etc.)
    'TEST_seed',47,...                      % seed value to for random number generator
    'TEST_target_rms',0.1,...                      % rms level for test sound
    'TEST_bandwidth',[4000 8000],...        % frequency range (min max)
    'TEST_dur',0.1,...                      % duration of test sounds (seconds)
    'TEST_ramp',5,...                       % ramp duration for test sound (should be 5ms)
    'TEST_on_delay_pts',[],...              % delay before test sound onset (TBD)
    'TEST_trial_freq',10, ...                % test sound every x trials
    'TEST_trial_jitter',5, ...              % maximum jitter in isi_buf
    'TEST_SPLs',[-12 0 3 6 12],...          % target SPLs for test sounds
    'TEST_scales',[],...                    % scales (calculated using calibration data)
    'TEST_sound',[],...                     % test sound stored here
    'TEST_scale_sequence',[],...            % sequence of test sound scales
    'TEST_nlocs',[],...                     % no. of test locations (determined from TEST_locs list)
    'TEST_locs',[0 -85; 0 -70; 0 -55; 0 -40; ...
    0 -25; 0 -10; 0 5; 0 0; 0 5; 0 10; 0 25; 0 40; 0 55; 0 70; 0 85],...      % locations for test stimuli (el,az)
    'TEST_loc_sequence',[],...              % sequence of test sound locations
    ...
    ...                                     % HRTF PARAMETERS:
    'HRTF_directory',hrtf_path,...  % directory of HRTF coefficient files
    'HRTF_fname',[],...
    'HRTF_nTaps',255,...
    ...
    'CALIB_directory',calib_path,...    % directory for calibration files
    'CALIB_fname',[]);

% calculate # locs
PDR.TEST_nlocs=size(PDR.TEST_locs,1);

% Select HRTF and Calib filenames:
switch(PDR.bird_id)
    case 1073
        PDR.HRTF_fname='1073AC_eq_ABLequal_normed.mat';
        PDR.CALIB_fname='...';
    case 925
        PDR.HRTF_fname='925AD_eq_ABLequal_normed.mat';
        PDR.CALIB_fname='intraural_calib_201358_925S.mat';
end

