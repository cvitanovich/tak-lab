% NAME: run_spatial_adapt.m
%
% AUTHOR: A Cvitanovich (incorporates some code by K Keller, L Whitchurch, K Hartung, and ADS Bala)
%
% CREATED: July 2012
%
% DESCRIPTION: Adapted from run_lizpdr_virt for a spatial adaptation pdr
% project. For presentation of a sequence of test probes with an adapting stimulus.
% Both the test probes and adapting stimulus can be presented from selected
% locations in acoustic space by convolving with individualized HRIRs.
% Freefield presentation has not been set up for the time being.
% This is the main script for sound output through the TDT system II:
%
% WIRING:
%
%
% DAC[0] -> Left PF1    -> Left PA4     -> Left In (HB6)    -> Left Earphone
% DAC[1] -> Right PF1   -> Right PA4    -> Right In (HB6)   -> Right Earphone
%
%
% NOTES: This is the main script for running the experiment. All dependent
% functions are kept in a "switchyard" called: spatialPDR_switchyard.m
% All important information about each experiment is kept in a global structure
% called simply "PDR" which is saved in the folder: C:\andrew\pdr\data\
% as header files which have names with the format YYMMDDxy_spatial.mat
% where x = bird's identification (l = #929, m = #930, etc.) and
% y = session ID (a, b, c, etc.).
%
% DEPENDENCIES:
%
% c:\andrew\pdr\code\spatialPDR.c -- Code for TDT control (adapted from play2_record2D2.c by K Keller).
% This needs to be compiled using mex.
%
% c:\andrew\play.c -- Functions used by spatialPDR.c, including the functions play0 and
% play0_record0 as well as FIR filtering functions.
%
% c:\andrew\play.h -- C header file.
%
% c:\andrew\pdr\code\spatialPDR_calib.c -- C code for calibration sessions.
%
% c:\andrew\play_calib.c -- Calibration functions;
%
% c:\andrew\play_calib.h -- Calibration header;
%
% c:\TDT\S232\VC\s2drv32c.lib -- C library for TDT control.
%
% spatialPDR_switchyard.m -- the function switchyard
%
% The following functions are in spatialPDR_switchyard.m:
%
% setDefaults() -- Sets up defaults for experiment!
%
% quit = setupTrialSeq() -- Sets up the trial sequence for each exp't!
% Outputs a "quit" flag if user decides to quit.
%
% MTLreadHDR() -- Reads HRTF file headers (original: MTLRH.m by K Hartung 1995)
%
% MTLreadDIR() -- Reads direction matrix from HRTF files 
% (original: MTLRDIR.m by K Hartung 1995)
%
% [dir1, dir2] = sphere2double(p1, p2) -- Converts spherical to double polor
% coordinates (original: SPH2DBL.m by K Hartung 1995). p1 and p2 are the
% spherical coordinates (takes either 1 or 2 arguments)
% and dir1/dir2 are the double polar outputs.
%
% [channel] = MTLreadCH(index) -- Reads channels from HRTF files
% (original: MTLRCH.m by K Hartung 1995). Input "index" is the index of the
% channel, and output "channel" is the desired channel.
%
% setupCalib(flag_adapt) -- Sets up a calibration session using longer
% sounds. Requires the S232 dll because the AP2 stack cannot handle the many
% large buffer sizes needed. We get around this by doing most computations
% in Matlab and sending APOS and XBUS commands using S232.dll. This permits
% longer sounds to be played, but the timing of each stimulus is not as
% precise as in experimental sessions. No matter, because this is not a concern for
% calibrations.
%
% stim = makeTest(seedval) -- Makes test sounds using specific parameters,
% including bandwidth, seed value, and RMS level. Input "seedval" is the
% seed value and output "stim" is the test sound produced.
%
%
% [X,Xtime] = seeded_whnoise(minfreq,maxfreq,Fs,duration,seedval,phi) --
% Just like whnoise.m but uses a specified seed value for reproducible sound.
%
% run_calibration(passmein) -- Runs a calibration session. "passmein" is a
% structure with session parameters.

button = questdlg('Clear all variables, close all windows and start?','Spatial Adaptation','YES','NO','YES');


if strcmp(button,'NO')
    return
else
    clear all; close all;
end

global PDR HRTF session

% Select to do experiment or calibrations:
type = menu('Experiment or calibrations?','Experiment  ','Calibrations');
if type == 2
    flg = menu('What do you wish to calibrate?','Test Probes','Adaptor');
    flag_adapt = flg-1;
    spatialPDR_switchyard('setupCalib',flag_adapt);
   
    % Setup HRTFs for calibration session:
    spatialPDR_switchyard('MTLreadHDR');
    spatialPDR_switchyard('MTLreadDIR');
    
    %convert coordinates to double polar, only if not using 929 or 930 or ones
    if (~strcmp(PDR.HRTF_fname(1:3),'930') & ~strcmp(PDR.HRTF_fname(1:3),'929'))
        PDR.HRTF_dir_matrix = spatialPDR_switchyard('sphere2double',PDR.HRTF_dir_matrix);
    end
    
    HRTF.CALIB_Left = nan*ones(PDR.HRTF_nlines,PDR.CALIB_nlocs);
    HRTF.CALIB_Right = HRTF.CALIB_Left;
    direc = PDR.HRTF_dir_matrix; % NOTE: first row = Elevation and 2nd row = Azimuth !!!
    for i=1:PDR.CALIB_nlocs
        idx{i}=find(direc(1,:)==PDR.CALIB_locations(1,i) & direc(2,:)==PDR.CALIB_locations(2,i));
        HRTF.CALIB_Left(:,i) = single(spatialPDR_switchyard('MTLreadCH',idx{i}*2-1));
        HRTF.CALIB_Right(:,i) = single(spatialPDR_switchyard('MTLreadCH',idx{i}*2));
    end

    % Run the calibration session:
    quit_flag = spatialPDR_switchyard('run_calibration');
    if quit_flag
        warndlg('Session Aborted!','Aborted Session');
        if exist([PDR.data_path PDR.filename '.mat'],'file')
            answer = questdlg('Delete session file?','Aborted Session','YES','NO','YES');
            if strcmp(answer,'YES')
                delete([PDR.data_path PDR.filename '.mat']);
                txt = ['Deleted file: ' PDR.data_path PDR.filename '.mat'];
                warndlg(txt,'File Deleted');
            end
        end
        close all;
    end
    return;
end

spatialPDR_switchyard('setDefaults'); % sets default values for running a session

quit = 0;
quit = spatialPDR_switchyard('setupTrialSeq'); % just sets up trial IDs so far

if quit
    return;
end

% RAMP (For Test Sound Buffer) SETUP:
ramplen = 5; % length of ramp in ms
tmp = find(PDR.TEST_sound ~= 0);
start = tmp(1);
stop = tmp(end);
ptsramp=round(ramplen/1000*PDR.stim_Fs);
on_rmp=1:-1/(ptsramp-1):0;
off_rmp=0:1/(ptsramp-1):1;
PDR.ADAPT_ramp = ones(1,PDR.buf_pts);
PDR.ADAPT_ramp(1,start-length(on_rmp):start-1)=on_rmp;
PDR.ADAPT_ramp(start:stop) = zeros(1,stop-start+1);
PDR.ADAPT_ramp(1,stop+1:stop+length(off_rmp))=off_rmp;

if 0 % not using a menu right now
H = spatialMenu; % setup menu
uiwait(H)

if PDR.exit_flag == -1
    disp('Starting experiment...')
elseif PDR.exit_flag == 1
    disp('Goodbye!')
    return
else
    warndlg('???')
    return
end
end

% HRTF SETUP:

%read HRTF coefficients files:
spatialPDR_switchyard('MTLreadHDR');
spatialPDR_switchyard('MTLreadDIR');

%convert coordinates to double polar, only if not using 929 or 930 or ones
if (~strcmp(PDR.HRTF_fname(1:3),'930') & ~strcmp(PDR.HRTF_fname(1:3),'929'))
    PDR.HRTF_dir_matrix = spatialPDR_switchyard('sphere2double',PDR.HRTF_dir_matrix);
end
HRTF.TestL = nan*ones(PDR.HRTF_nlines,PDR.TEST_nlocs);
HRTF.TestR = HRTF.TestL;
direc = PDR.HRTF_dir_matrix; % NOTE: first row = Elevation and 2nd row = Azimuth !!!
for i=1:PDR.TEST_nlocs
    idx{i}=find(direc(1,:)==PDR.TEST_locs(i,1) & direc(2,:)==PDR.TEST_locs(i,2));
    HRTF.TestL(:,i) = single(spatialPDR_switchyard('MTLreadCH',idx{i}*2-1));
    HRTF.TestR(:,i) = single(spatialPDR_switchyard('MTLreadCH',idx{i}*2));
end

if PDR.flag_adapt
    % Make FIR coefficients for Adaptor:
    HRTF.AdaptL = nan*ones(1,PDR.HRTF_nlines);
    HRTF.AdaptR = nan*ones(1,PDR.HRTF_nlines);
    spatialPDR_switchyard('makeGammaFIR');  
    idxm=find(direc(1,:)==PDR.ADAPT_loc(1,1) & direc(2,:)==PDR.ADAPT_loc(1,2));
    % HRTFs for adapting stimulus:
    HRTF.AdaptL(1,:) = single(spatialPDR_switchyard('MTLreadCH',idxm*2-1));
    HRTF.AdaptR(1,:) = single(spatialPDR_switchyard('MTLreadCH',idxm*2));
else
    HRTF.AdaptL = zeros(1,PDR.HRTF_nlines);
    HRTF.AdaptR = HRTF.AdaptL;
end

% Make sure the session gets a unique filename
cd(PDR.data_path);
cnt = double('a'+0);
while exist ([PDR.filename '.mat'],'file');
    cnt = cnt + 1;
    if cnt > 122
        disp(['Attempted to use filename: ' PDR.filename '.mat' ' but failed!']);
        disp('There are already several files with similar names!');
        PDR.filename = input('Enter a unique filename for this session: ', 's');
        break;
    else
        PDR.filename(end) = char(cnt);
    end
end

button1 = questdlg(['Check that equipment is configured correctly.' ...
        'Is it okay to continue with this session?'],'','YES','NO','YES');

if strcmp(button1,'NO')
    return
end



%*********************************%
%** PLOT PDR TRACE ***************%
%*********************************%
session.bufpts = PDR.buf_pts;
session.dec_fact = PDR.decimationfactor;
session.isi = PDR.isi_buf;
session.trials_to_show = 3;
tmp=find(PDR.TEST_sound~=0);
session.sound_onset=tmp(1); % how many points until sound onset?
session.npts_totalplay=PDR.npts_totalplay;
session.srate=(10^6)*(1/PDR.stim_Fs); % sampling period in usec
session.stim_fs=PDR.stim_Fs;
session.zoomval=0.4;
screen_size = get(0, 'ScreenSize');% get scrn size
session.hTracePlot=figure; whitebg(gcf,'k');
set(session.hTracePlot,'renderer','OpenGL'); %use OpenGL for renderer
set(session.hTracePlot, 'Position', [0.5*screen_size(3) 0.04*screen_size(4) 0.45*screen_size(3) 0.45*screen_size(4)] );
hold on;

% requires session to be declared a global
cd ..\code
sessionPlots2('Initialize');
hold on;
uicontrol(session.hTracePlot,'Style', 'pushbutton','Tag','ZoomOut','String','Zoom -',...
    'Units','normalized','FontSize',8,'Position',[0.02 0.6 0.05 0.05],...
    'Callback', 'if session.zoomval<2.1; session.zoomval=session.zoomval+0.1; end;');
uicontrol(session.hTracePlot,'Style', 'pushbutton','Tag','ZoomIn','String','Zoom +',...
    'Units','normalized','FontSize',8,'Position',[0.02 0.8 0.05 0.05],...
    'Callback', 'if session.zoomval>0.1; session.zoomval=session.zoomval-0.1; end;');
% required to initialize sessionPlots:
% structure with these parameters:
% hab_xes, hab_yes, trial_xes, trial_yes
% ntrials, min_yes, max_yes
% bufpts, decpts, isi, trials_to_show

%*********************************************************************
%*********************************************************************
%       Start Main Loop Here
%*********************************************************************
%*********************************************************************
%from here on, this will be run in a .dll and you can't jump out.

lengthOFtrials=(round(PDR.npts_totalplay/PDR.stim_Fs/60*10))/10;  %should be in approximate minutes
b=clock;
PDR.starttime(1:2)=b(4:5);
disp(['Now the program will just have to run its course.'])
disp(['It should be done in ~' num2str(lengthOFtrials) ' min from now (' num2str(b(4:5)) ')']);

passmein.outFN1=[PDR.filename '_REC1.vrt'];
passmein.outFN2=[PDR.filename '_REC2.vrt'];
passmein.record = PDR.record;
passmein.n_trials = PDR.ntrials;
passmein.latten = PDR.base_atten;
passmein.ratten = PDR.base_atten;
passmein.adaptScale = PDR.ADAPT_scale;
passmein.decimateFactor = PDR.decimationfactor;
passmein.buf_pts = PDR.buf_pts;
passmein.nptsTotalPlay = PDR.npts_totalplay;
passmein.ISI = PDR.isi_buf;
passmein.flag_adapt = PDR.flag_adapt;
passmein.max_signal_jitter = PDR.TEST_trial_jitter;
passmein.n_locs = PDR.TEST_nlocs;
passmein.HRTF_nlines = PDR.HRTF_nlines;

cd(PDR.data_path)
p = pwd;
p = [p '\'];
if strcmp(PDR.data_path,p)
else
    warndlg('Something could be wrong with the path setup!')
    return
end

%write header information to file... saving global variables
PDR = orderfields(PDR); % order fields by ASCII dictionary order
save ([PDR.data_path PDR.filename '.mat'], 'passmein','PDR','HRTF');
%RUN dll to control TDT
spatialPDR(passmein, single(HRTF.AdaptL),single(HRTF.AdaptR), single(HRTF.TestL),single(HRTF.TestR), uint32(PDR.ADAPT_seed),single(PDR.ADAPT_coefs), single(PDR.ADAPT_ramp), single(PDR.TEST_sound), single(PDR.TEST_scale_sequence), int32(PDR.TEST_loc_sequence));

b=clock;
PDR.stoptime(1:2)=b(4:5);
%write header information to file... saving global variables
cd(PDR.data_path)
save ([PDR.data_path PDR.filename '.mat'], 'passmein','PDR','HRTF');
cd ..\code