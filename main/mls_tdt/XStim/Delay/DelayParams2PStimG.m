function [PStimG] = DelayParams2PStimG(XStimParams, FN);

PStimG = struct('TestName', 'Delay',...
    'TestNameNum', XStimParams.testnum,...
    'Note', 'XStim',...                          % new - note
    'DataPath', FN.data_path,...                
    'RCOpath', 'C:\TDT\Brian\',...              % new - path to RCO files
    'RCOfile', 'SysII',...                      % new name of DEFAULT .rco file
    'RCOfileNum', NaN,...                       % new - file number from list
    'GlobalPath', FN.home,...                   % path to globals file
    'Sys3Interface', 'SysII',...                % Device interface type
    'DevNum', 1,...                             % new - Device number
    'isConnected', 0,...                        % new - Processor connected correctly
    'SR', 30000,...                             % TDT variable not saved - Device sampling rate
    'HRTFfile', FN.space,...
    'HRTFpath', FN.space_path,...
    'HRTFformat', FN.HRTFfiletype,...           % ???? HRTF file format
    'HRTFdelay', 1.5,...                        % new - HRTF delay estimate
    'SaveData', 1,...                           % new - not saved but obviously yes, Save data
    'SpaceResolution', XStimParams.GUIspaceres,...   % Space picker resolution
    'SimpleRaster', 1,...                       % new - plot simple raster
    'CursorMode', XStimParams.picker_flag,...        % ???  Space picker cursor mode
    'Session', XStimParams.session_num,...           % Session string
    'BirdNum', XStimParams.bird_number,...           % Bird number
    'RecordingSite', XStimParams.recording_site,...  % Recording site
    'TestNum', XStimParams.testnum,...               % Test number
    'AP', XStimParams.ap_pos,...                     % AP position
    'ML', XStimParams.ml_pos,...                     % ML position
    'Depth', XStimParams.depth,...                   % Depth
    'Locations', XStimParams.locations,...           % El and Az locations
    'Az_1', XStimParams.locations(1),...             % Best location
    'El_1', XStimParams.locations(2),...             % Best location
    'Az_2', XStimParams.offset_az,...                % Second source
    'El_2', XStimParams.offset_el,...                % Second source
    'DelayMods', XStimParams.DelayMods,...           % Delay modifications (see other code)
    'StimDelays', XStimParams.DelayTimes,...         % Delay delays
    'OnsetRamps', XStimParams.ramp_timeS_on,...      % Delay onset ramp times
    'OffsetRamps', XStimParams.ramp_timeS_off,...    % Delay offset ramp times
    'RampStim', XStimParams.OnOff_mode,...           % Delay ramp lag|lead & lag|lead
    'LagExts', XStimParams.lag_seg_ext,...           % Delay lag segment extensions
    'Atten', XStimParams.curr_ABL,...                % Attenuation
    'HBx_Atten', 0,...                          % new - effectively zero
    'Atten_Steps', [XStimParams.loabl (XStimParams.hiabl-XStimParams.loabl)/XStimParams.numabls XStimParams.numabls],... % min attenuation (dB), step, # steps for Rate Level
    'Freq_Steps', [XStimParams.lofreq (XStimParams.hifreq-XStimParams.lofreq)/XStimParams.numfreqs XStimParams.numfreqs],...  % min frequency (kHz), step, # steps for Rate Freq
    'StimDur', XStimParams.curr_stimdur,...          % Stimulus Duration
    'PreStimDelay', XStimParams.silence_lead,...     % time spikes collected before stim
    'PostStimDelay', XStimParams.silence_trail,...   % time spikes collected after stim
    'ISI', XStimParams.test_ISI,...                  % Interstimulus Interval
    'NumReps', XStimParams.numreps,...               % Number repetitions
    'ITD', XStimParams.curr_ITD,...                  % ITD
    'ILD', XStimParams.curr_ILD,...                  % ILD
    'addReverb', NaN,...                          % new - add reverb (if reverb components are in chain)
    'WetMix', NaN,...                             % new -Reverb WetMix
    'Decay', NaN,...                            % new - Reverb Decay
    'Dcurr', NaN,...                             % new - Reverb Delay
    'UseHRTF', XStimParams.space_flag,...            % Convolve with HRTFs
    'PlayTone', NaN,...                           % new but indicated by XStimParams.stim_type1..3
    'Frequency', XStimParams.curr_freq,...           % frequency of tone
    'PlotStim', NaN,...                           % new - plot stimuli during trial
    'duh', NaN);                                  % new -end
