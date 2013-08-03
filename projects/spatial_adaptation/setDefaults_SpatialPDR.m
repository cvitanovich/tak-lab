% SETUP DEFAULTS FOR spatialPDR EXPERIMENT

% select folder location based on computer being used:
if(~ispc && strcmp(getenv('USER'),'cvitanovich'))
    code_path='/Users/cvitanovich/Documents/MATLAB/tak-lab/projects/spatial_adaptation/';
    data_path='/Users/cvitanovich/Documents/MATLAB/data/';
    hrtf_path='/Users/cvitanovich/Documents/MATLAB/tak-lab/HRTFs/Matlab_V6/';
    calib_path='/Users/cvitanovich/Documents/MATLAB/tak-lab/calib_data/spatial_PDR/';
else
    code_path='C:\Documents and Settings\andrew\My Documents\GitHub\tak-lab\projects\spatial_adaptation\';
    data_path='C:\andrew\data\';
    sounds_path='C:\andrew\sounds\';
    hrtf_path='C:\Documents and Settings\andrew\My Documents\GitHub\tak-lab\HRTFs\Matlab_V6\';
    calib_path='C:\Documents and Settings\andrew\My Documents\GitHub\tak-lab\calib_data\spatial_PDR\';
end
    
PDR = struct(...                            % MAIN PARAMETERS:
    'DEBUG',0,...                           % set to zero unless debugging without TDT!
    'jitter',0,...
    'virtual',1, ...                        % flag for virtual sound presentation
    'record',1, ...                        % flag for recording pupillometer output
    'flag_adapt',1,...                     % flag for adapting stimulus
    'ntrials',10,...                       % # of trials in session
    'npretrials',2,...                     % # of trials before first test trial
    'n_test_trials',[],...                  % # of test trials 
    'buf_pts',30000,...                     % number of pts in each buffer
    'buf_dur',[],...                        % buffer duration in ms
    'isi_buf',3,...                         % no. buffers between trials
    'isi_time',[],...                       % ISI (seconds) between trials
    'decimationfactor',4, ...               % decimation factor for data collection
    'stim_Fs',30000, ...                   % sampling rate (Hz) for stimuli
    'npts_totalplay',[],...                % total no. pts in the session
    'len_session',[], ...                   % length of session (in minutes)
    'starttime',[], ...                     % session start time
    'stoptime',[], ...                     % session stop time
    'code_path',code_path,...  % path to code
    'data_path',data_path, ...     % data path (normally: c:\andrew\pdr\data normally, calib: c:\andrew\pdr\calib_data)
    'base_atten',40,...70,...                      % base attenuation for PA4
    'filename',[],...                       % filename for session data
    'bird_id',925,...
    'hab_id',0,...
    'sounds_path',sounds_path,...           % path to sound files
    ...
    ...                                     % ADAPTING STIMULUS PARAMETERS:
    'ADAPT_file',[sounds_path 'gtone_B.mat'],... % file with adaptor sound
    'ADAPT_loc',[0 0],...                   % location of adapting stimulus
    'ADAPT_species','owl',...               % species
    'ADAPT_soundtype','GTONE',...            % type of sound for adapting stimulus
    'ADAPT_cF',6000,...                     % center frequency (Hz) (if relevant)
    'ADAPT_coefs',[],...                    % FIR coefficients to generate adapting stimulus
    'ADAPT_SPL',30,...                      % SPL for adaptor
    'ADAPT_scale',[],...                    % scale for adaptor (calculated using calibration data)
    'ADAPT_ramp',5,...                      % ramp for trial segments
    'ADAPT_dur',[],...                      % adaptor duration in seconds (calculated from isi_buf) in seconds
    'ADAPT_target_rms',0.1,...              % desired rms amplitude for adaptor (before scaling)
    'ADAPT_seeds',[],...                    % seed for generating adaptor
    'ADAPT_sound',[],...
    ...
    ...                                     % TEST STIMULUS PARAMETERS:
    'TEST_soundtype','OCTAVE',...           % type of test sound ('gammatone', 'octave', 'broadband', etc.)
    'TEST_seed',47,...                      % seed value to for random number generator
    'TEST_target_rms',0.1,...                      % rms level for test sound
    'TEST_bandwidth',[4000 8000],...        % frequency range (min max)
    'TEST_dur',0.1,...                      % duration of test sounds (seconds)
    'TEST_ramp',5,...                       % ramp duration for test sound (should be 5ms)
    'TEST_on_delay_pts',[],...              % delay before test sound onset (TBD)
    'TEST_trial_freq',3, ...                % test sound every x trials
    'TEST_trial_jitter',0, ...              % maximum jitter in isi_buf
    'TEST_SPLs',[20 30 40],...[-12 -6 -3 0 3 6 12],...          % target SPLs for test sounds
    'TEST_scales',[],...                    % scales (calculated using calibration data)
    'TEST_sound',[],...                     % test sound stored here
    'TEST_scale_sequence',[],...            % sequence of test sound scales
    'TEST_loc',[0 0],...                   % location for test stimuli (el,az)
    ...
    ...                                     % HRTF PARAMETERS:
    'HRTF_directory',hrtf_path,...  % directory of HRTF coefficient files
    'HRTF_fname',[],...
    'HRTF_nTaps',255,...
    ...
    'CALIB_directory',calib_path,...    % directory for calibration files
    'CALIB_fname',[]);

% Select HRTF and Calib filenames:
switch(PDR.bird_id)
    case 1073
        PDR.HRTF_fname='1073AC_eq_ABLequal_normed.mat';
        PDR.CALIB_fname='...';
    case 925
        PDR.HRTF_fname='925AD_eq_ABLequal_normed.mat';
        PDR.CALIB_fname='intraural_calib_201358_925S.mat';
end

