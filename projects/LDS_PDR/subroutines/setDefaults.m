function setDefaults()
global PDR

if(~ispc)
    code_path='/Users/cvitanovich/Documents/MATLAB/tak-lab/projects/LDS_PDR/';
    data_path='/Users/cvitanovich/Documents/MATLAB/data/';
    calib_path='/Users/cvitanovich/Documents/MATLAB/calib/';
else
    code_path='c:\alex\code\LDS_PDR\';
    data_path='C:\alex\data\';
    calib_path='C:\alex\calib_data\';
end
cd(code_path);

PDR = struct(...                        % MAIN PARAMETERS:
    'DEBUG',1,...
    'bird_id',1073,...                  % bird's id #
    'info_pts',2,...                    % first 2 points in the (decimated) buffer reserved for trial info
    ...                                 % for example, if testing MAAs a comment should go here
    'maa_test_flag',0,...               % set to one if testing MAAs            
    'virtual',0,...                     % Freefield: virtual = 0, Headphones: virtual = 1
    'record',1,...                      % flag for recording pdr trace
    'jitter',2,...                      % jitter trial occurance (2 bufs = +/- 1.34 s)
    'ntrials',[],...                   % no. of trials in a session
    'npretrials',100,...                 % default # of habituating trials before the 1st test trial
    'n_test_trials',15,...              % default # of test trials in a session
    'buf_pts',32768,...                 % # points in a buffer - 32768 pts for a 0.671 second buffer (Fs = 48828 Hz)
    'buf_dur',[],...                    % buffer duration in ms
    'isi_buf',15,...                     % number of empty buffers between trials (10 buffers gives an ISI of ~6.71s (+/-1.34 s) using .671 sec sounds)
    'isi_time',[],...                   % time between trial buffers in seconds (calculated)
    'decimationfactor',5,...            % decimation factor for data collection (using 4 will give 1875Hz sampling rate)
    'stim_Fs',48828,...                 % sampling rate for sound production
    'npts_totalplay',[],...            % total no. pts in the session
    'len_session',[nan nan], ...               % length of session (in minutes)
    'starttime',[], ...                 % session start time
    'stoptime', [], ...                 % session stop time
    'code_path',code_path,...    % path to code
    'data_path',data_path,...     % flag indicates that AD recording (pupillometer) will be needed
    'base_atten',0,...                  % attenuation value to send to TDT
    'filename',[],...                   % file name for storing session data
    'exit_flag',0,...                   % exit session (quit) flag
    'plot_trials_handle',[],...                 % handles for the menu plots
    'plot_sounds_handle',[],...
    'outlier_separation',100,...        % defines angular separation considered to be "extreme outlier"
    'outlier_fraction',[],...          % fraction of test trials that are outliers (10% seems fine)
    'n_outlier_trials',0,...            % fraction of tests that are outliers
    'RPs',[],...                        % store response probabilities (From Nelson AM study)
    ...
    ...                                 % LEAD SOUND PARAMETERS:
    'LEAD_pos',-60,...                  % default azimuth (degrees) for lead sound
    'LEAD_sounds',[],...                 % LEAD_sounds sound (to be calculated)
    'LEAD_default_scale',0,...          % will be set to zero for MAA testing
    ...                                 % otherwise, lead/lag scales should be identical
    ...
    ...                                 % LAG SOUND PARAMETERS:
    'LAG_hab_pos',-50,...                % default azimuth for lag sound
    'LAG_sounds',[],...                  % lag sound (to be calculated)
    'LAG_default_scale',0,...       % lead and lag should have same scale (unless testing MAAs)
    ...
    ...                                 % TEST SOUND PARAMETERS:
    'TEST_azimuths',[-25 -5 5 15 30],... %45 60],... % SHOULD BE IN ORDER FROM SMALLEST TO GREATEST SEPARATION!!!!
    'TEST_trial_jitter',5,...           % jitter (test trials)
    'TEST_trial_freq',30,...             % default frequency of test trials (1 out of every X trials)
    ...
    ...                                 % SOUND PARAMETERS:
    'CARRIERS_SWITCHED',0,...           % if the delay is negative this flag will be set to one! (LEAD <-> LAG carriers switched if one)
    'SOUNDS_calib_path',calib_path,...  % path to files used for setting calibrated attens/scales
    'SOUNDS_rms',0.5,...                % setting envelopes to 0.5 rms amplitude (equalizes sound envelope magnitudes)
    'SOUNDS_radius',152,...                % the distance from the speakers to the owl's head position (in cm)
    'SOUNDS_calib_fname',[],...         % filename for speaker calibration parameters
    'SOUNDS_SPL',60,...                  % target decibel level (intraural) for speaker
    'SOUNDS_num_speakers',7,...         % number of speakers (including lead speaker) available
    'SOUNDS_speaker_scales_lead',[],...      % list of calibrated scale values for lead speaker (matched to each lag)
    'SOUNDS_speaker_scales_lag',[],...      % list of calibrated scale values for lag speaker
    'SOUNDS_lead_attens',[],...             %attens for lead speaker
    'SOUNDS_elevations',[0 0 0 0 0 0 0],... %0 0],...
    'SOUNDS_speaker_numbers',1:7,...
    'SOUNDS_azimuths',[-60 -50 -25 -5 5 15 30],... % 45 60],...            % actual speaker azimuths (don't change unless moving speakers)
    'SOUNDS_length',100,...              % length of sounds in ms
    'SOUNDS_carrier_delay',-3,...                % ongoing CARRIER disparity between lead/LAG_sound sounds
    'SOUNDS_env_delay',0,...            % ongoing ENVELOPE disparity (lead/lag)
    'SOUNDS_mix',0,...                % mix fraction
    'SOUNDS_env_correlations',[],...     % calculated correlation coefficient (envelopes)
    'SOUNDS_states',[795 4932 4584],...% state values for random no generator (to ensure that sounds are reproducible)
    ...                                 % formatted as such: [env_state decorr_env_state1 decorr_env_state2]
    'SOUNDS_carrier_states',[],...            % carrier states (to rove the carrier)
    'SOUNDS_num_carriers',10,...        % number of carriers to rove (should be 10, don't change!)
    'SOUNDS_rand_states',[],...         % list of random state values (from Caitlin/Brian head saccade exp't)
    'SOUNDS_rove_sequence',[],...      % random ordering sequence to rove sounds
    'SOUNDS_carrier_bandwidth',[2000 11000],... % using 2 - 11 kHz bandwidth for carriers (same as Caitlin/Brian head saccade exp't)
    'SOUNDS_env_bandwidth',[1 150],...          % bandwidth of envelopes (should be 1 to 150 Hz to replicate Nelson/Baxter study)
    'SOUNDS_amplitude',60,...                   % Keep at 60 (this gives amplitudes suitable for the TDT)
    'SOUNDS_env_depth',100,...                  % depth of envelope modulation
    'SOUNDS_ramp',5,...                         % ramp length for sounds (ms)
    'SOUNDS_location_sequence',[]);           % sequence of locations for the session (1st row: lead, 2nd row: lag)

% Comments box:
PDR.comments = cell(3,1);
for j=1:length(PDR.comments)
    PDR.comments{j}='';
end
% These are the random SOUNDS_states from Baxter, et al.
RandStates = [...
        73 167 617 364 662 593 538 194 853 610 294 ...
        479 71 105 162 770 143 116 252 101 377 706 ...
        273 574 915 661 935 355 530 540 220 232 886 ...
        70 65 571 35 339 87 281 795 283 974 248 ...
        995 936 769 943 127 224];   
PDR.SOUNDS_rand_states=RandStates;

% check if testing MAAs
if PDR.maa_test_flag
    PDR.LEAD_default_scale = 0;
end

% setup cell array for lead/lag sounds
PDR.LEAD_sounds=cell(1,PDR.SOUNDS_num_carriers);
PDR.LAG_sounds=cell(1,PDR.SOUNDS_num_carriers);

q=clock;
y=num2str(q(1));y=y(3:4);
m=num2str(q(2));if size(m,2)<2;m=['0' m];end
d=num2str(q(3));if size(d,2)<2;d=['0' d];end
PDR.filename= [num2str(PDR.bird_id) '_LDS_' y m d '_a'];                    %930=m 929=l 882=d 883=e  (a is index to experiment number for that day)