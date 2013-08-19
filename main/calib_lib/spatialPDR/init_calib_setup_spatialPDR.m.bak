function init_calib_setup_spatialPDR
% initialize TDT and set up for calibrations
global PDR H TDT KNOWLES

% CHECK IF IT'S OKAY TO CONTINUE
button = questdlg('You need to reset the TDT before continuing. Click OKAY when you''ve done this, or CANCEL otherwise.',...
'Continue Operation','OKAY','CANCEL','OKAY');
if strcmp(button,'OKAY')
elseif strcmp(button,'CANCEL')
   return;
end

% SETUP PARAMETERS
% NOTE: This structure should be saved as CALIB_PDR to avoid confusion
% but for all the scripts to work it needs to be called PDR for now!

% select folder location based on computer being used:
if(ismac && strcmp(getenv('USER'),'cvitanovich'))
    code_path='/Users/cvitanovich/Documents/MATLAB/tak-lab/projects/calib_lib/spatialPDR/';
    data_path='/Users/cvitanovich/Documents/MATLAB/data/calib/';
    hrtf_path='/Users/cvitanovich/Documents/MATLAB/tak-lab/HRTFs/Matlab_V6/';
else
    code_path='C:\andrew\CORE\tak-lab\calib_lib\spatialPDR\';
    data_path='C:\andrew\calib_data\';
    hrtf_path='C:\andrew\CORE\tak-lab\HRTFs\Matlab_V6\';
end

PDR = struct(...                        % MAIN PARAMETERS:
    'calib',1,...                       % flag for calibrations
    'owl_id',0,...                      % for owl ID #
    'virtual',1,...                     % Freefield: virtual = 0, Headphones: virtual = 1
    'record',1,...                      % flag for recording CALIB_PDR trace
    'buf_pts',30000,...              
    'stim_Fs',30000,...                 % sampling rate for sound production (1 second calibration sounds)
    'npts_totalplay',[],...             % total no. pts in the session
    'len_session',[], ...               % length of session (in minutes)
    'starttime',[], ...                 % session start time
    'stoptime', [], ...                 % session stop time
    'code_path',code_path,...    % path to code
    'data_path',data_path,...
    'base_atten',0,...                  % attenuation value to send to TDT (combining scale_val and atten gives a certain SPL in dB -- This needs to be calibrated!!!)
    'filename',[],...                   % file name for storing session data
    ...
    ...                                 % CALIBRATION PARAMETERS:
    'TEST_ramp',5,...
    'TEST_seed',47,...
    'TEST_base_rms',0.0350,...
    'calibtype',[],...                  % What are we calibrating? (Knowles alone, Knowles in owl's ears)
    'knowles_file',[],...               % equals 1 if calibrating Knowles voltage output
    'comments','Spatial PDR with Octave Band, BBN and Gammatones (virtual). Intraural Calibrations with Knowles mics.',...             % extra info about sounds
    'nscales',300,...                   % number of scales used for calibrating each sound (5 minutes per sound token)
    'scales_2_try_for_cutoffs',fliplr(10.^([0:.09:4.5 log10(32760)])),... % these are the scales to try to avoid nonlinearities due to speaker limitations and noise floor, etc.
    ...
    'HRTF_directory',hrtf_path,...  % directory of HRTF coefficient files
    'HRTF_fname','1073AC_eq_ABLequal.mat');

Prompt='Enter OWL ID # for this calibration session';
Title='Owl ID';
Answer = inputdlg(Prompt,Title);
PDR.owl_id=str2num(Answer{1});
set(H.owl_id2,'String',Answer{1},'BackgroundColor','w');
c=clock;
PDR.filename = ['intraural_calib_' num2str(c(1)) num2str(c(2)) num2str(c(3)) '_' num2str(PDR.owl_id) 'A'];
count = double('A'+0);
while exist ([PDR.data_path PDR.filename '.mat'],'file');
    count = count + 1;
    if count > 90
        disp('There are already several files with similar names!');
        PDR.filename = input('Enter a unique filename for this session: ', 's');
        break;
    else
        PDR.filename(end) = char(count);
    end
end
set(H.intraural2,'String',[PDR.data_path PDR.filename '.mat']);


% MAKE CALIB SOUNDS
HRTF = make_calib_sounds;

%write header information to file... saving global variables
set(H.calib_snd_btn,'Enable','on');
set(H.owl_id2,'Enable','on','BackgroundColor','w','String',num2str(PDR.owl_id));
set(H.load_knowles_data_btn,'Enable','on');
CALIB_PDR=PDR;
save ([CALIB_PDR.data_path CALIB_PDR.filename '.mat'],'CALIB_PDR');
str{1} = 'Variables saved in: ';
str{2} = [CALIB_PDR.data_path CALIB_PDR.filename '.mat'];
clear CALIB_PDR
hMsg=msgbox(str); uiwait(hMsg);
