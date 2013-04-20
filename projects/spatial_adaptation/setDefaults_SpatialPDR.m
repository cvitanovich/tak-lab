% SETUP DEFAULTS FOR spatialPDR EXPERIMENT

PDR = struct(...                            % MAIN PARAMETERS:
    'calib',0,...                           % flag for calibrations
    'virtual',1, ...                        % flag for virtual sound presentation
    'record',1, ...                        % flag for recording pupillometer output
    'flag_adapt',1,...                     % flag for adapting stimulus
    'ntrials',50,...                       % # of trials in session
    'npretrials',5,...                     % # of trials before first test trial
    'n_test_trials',[],...                  % # of test trials 
    'buf_pts',16384,...                     % number of pts in each buffer
    'buf_dur',[],...                        % buffer duration in ms
    'isi_buf',3,...                         % no. buffers between trials
    'isi_time',[],...                       % ISI (seconds) between trials
    'decimationfactor',4, ...               % decimation factor for data collection
    'stim_Fs', 30000, ...                   % sampling rate (Hz) for stimuli
    'npts_totalplay',[],...                % total no. pts in the session
    'len_session',[], ...                   % length of session (in minutes)
    'starttime',[], ...                     % session start time
    'stoptime', [], ...                     % session stop time
    'code_path', 'C:\andrew\CORE\tak-lab\projects\spatial_adaptation\',...  % path to code
    'data_path', 'C:\andrew\pdr\data\', ...     % data path (normally: c:\andrew\pdr\data normally, calib: c:\andrew\pdr\calib_data)
    'base_atten',0,...                      % base attenuation for PA4
    'filename',[],...                       % filename for session data
    'bird_id',929,...
    ...
    ...                                     % ADAPTING STIMULUS PARAMETERS:
    'ADAPT_loc',[0 0],...                   % location of adapting stimulus
    'ADAPT_dB',[],...                       % TBD
    'ADAPT_species','owl',...               % species
    'ADAPT_type','gammatone',...            % type of sound for adapting stimulus
    'ADAPT_cF',6000,...                     % center frequency (Hz) (if relevant)
    'ADAPT_coefs',[],...                    % FIR coefficients to generate adapting stimulus
    'ADAPT_scale',30000,...                  % test scale for adaptor
    'ADAPT_ramp',[],...                     % ramp for trial segments
    'ADAPT_seed',47,...                     % seed for generating adaptor (using C code)
    ...
    ...                                     % TEST STIMULUS PARAMETERS:
    'TEST_soundtype','octave',...           % type of test sound ('gammatone', 'octave', 'broadband', etc.)
    'TEST_seed',47,...                      % seed value to for random number generator
    'TEST_base_rms',0.0350,...              % rms level for test sound
    'TEST_bandwidth',[4000 8000],...            % frequency range (min max)
    'TEST_dur',100,...                      % duration of test sounds (ms)
    'TEST_ramp',5,...                       % ramp duration for test sound (should be 5ms)
    'TEST_on_delay_pts',[],...              % delay before test sound onset (TBD)
    'TEST_trial_freq',3, ...                % test sound every x trials
    'TEST_trial_jitter',0, ...              % maximum jitter in isi_buf
    'TEST_scales',[1000 5000 15000 32760],...                 %[10 20 60 80 500],... % test scales (try 10, 20 ,40, 60, 80, 200, 300, 500, 1000)
    'TEST_sound',[],...                     % test sound stored here
    'TEST_scale_sequence',[],...            % sequence of test sound scales
    'TEST_nlocs',5,...                     % no. of test locations (default is 3)
    'TEST_locs',[0 0; 0 20; 0 40; 0 60; 0 80],...      % locations for test stimuli (el,az)
    'TEST_loc_sequence',[],...              % sequence of test sound locations
    ...
    ...                                     % HRTF PARAMETERS:
    'HRTF_directory','C:\andrew\CORE\tak-lab\HRTFs\Matlab_V6\',...  % directory of HRTF coefficient files
    'HRTF_fname','1073AC_eq_ABLequal.mat');

% buffer duration (ms):
PDR.buf_dur = 1000*PDR.buf_pts/PDR.stim_Fs;

% Time (in seconds) between test stimuli:
PDR.isi_time = ( PDR.buf_dur*PDR.isi_buf+(PDR.buf_dur-PDR.TEST_dur)/2 )/1000;

q=clock;
y=num2str(q(1));y=y(3:4);
m=num2str(q(2));if size(m,2)<2;m=['0' m];end
d=num2str(q(3));if size(d,2)<2;d=['0' d];end
switch PDR.bird_id
    case 929
        LETTER = 'l';
    case 930
        LETTER = 'm';
    case 0
        LETTER = 't';
end

PDR.filename = [y m d LETTER 'a'];  %930=m 929=l 882=d 883=e  (a is index to experiment number for that day)
PDR.npts_totalplay = PDR.ntrials*(PDR.isi_buf+1)*PDR.buf_pts; % Calculate length of session!
PDR.len_session = (1/60)*(PDR.npts_totalplay/PDR.stim_Fs); % length of session in minutes
h=msgbox(['Session will last approximately ' num2str(PDR.len_session) ' minutes']);
uiwait(h)
cd(PDR.code_path);

% make test sound:
stim = makeTest(PDR.TEST_seed,PDR.TEST_dur,PDR.TEST_bandwidth(1),PDR.TEST_bandwidth(2),PDR.stim_Fs,PDR.TEST_ramp,PDR.TEST_base_rms);
PDR.TEST_sound = zeros(1,PDR.buf_pts);
on_delay_pts = floor((PDR.buf_pts - length(stim))/2);
PDR.TEST_on_delay_pts = on_delay_pts;
PDR.TEST_sound(on_delay_pts+1:on_delay_pts+length(stim)) = stim; % place stimulus in buffer mid-section
PDR.TEST_sound = (0.999)*PDR.TEST_sound ./ (max(abs(PDR.TEST_sound)));
clear stim;