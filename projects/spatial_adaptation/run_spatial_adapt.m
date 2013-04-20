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
% DEPENDENCIES:
%
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
% stim = makeTest(seedval) -- Makes test sounds using specific parameters,
% including bandwidth, seed value, and RMS level. Input "seedval" is the
% seed value and output "stim" is the test sound produced.
%
%
% [X,Xtime] = seeded_whnoise(minfreq,maxfreq,Fs,duration,seedval,phi) --
% Just like whnoise.m but uses a specified seed value for reproducible sound.

button = questdlg('Clear all variables, close all windows and start?','Spatial Adaptation','YES','NO','YES');


if strcmp(button,'NO')
    return
else
    clear all; close all;
end

global PDR HRTF TDT session

setDefaults; % sets default values for running a session

quit = 0;
quit = setupTrialSeq; % just sets up trial IDs so far

if quit
    return;
end

% RAMP SETUP:
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

if PDR.exit_flag == -1
    disp('Starting experiment...')
elseif PDR.exit_flag == 1
    disp('Goodbye!')
    return
else
    warndlg('???')
    return
end

% HRTF SETUP:

HRTF.TestL = nan*ones(PDR.HRTF_nlines,PDR.TEST_nlocs);
HRTF.TestR = HRTF.TestL;
if PDR.flag_adapt
    HRTF.AdaptL = nan*ones(1,PDR.HRTF_nlines);
    HRTF.AdaptR = nan*ones(1,PDR.HRTF_nlines);
else
    HRTF.AdaptL = zeros(1,PDR.HRTF_nlines);
    HRTF.AdaptR = HRTF.AdaptL;
end

if(strcmp(PDR.HRTF_fname((end-3):end)=='.mat')) % DOT MAT FORMAT
    LT=zeros(1,PDR.HRTF_nlines);
    RT=LT;
    for i=1:PDR.TEST_nlocs
        [LT, RT] = readHRTFdotMAT(fname,EL,AZ);
        HRTF.TestL(:,i)=LT';
        HRTF.TestR(:,i)=RT';
    end
    if PDR.flag_adapt
        [HRTF.AdaptL(1,:), HRTF.AdaptR(1,:)] = readHRTFdotMAT(fname,EL,AZ);
    end
else % using another format (e.g. for 930 or 929)
    %read HRTF coefficients files:
    [HRTF] = MTLreadHDR(PDR.HRTF_directory,PDR.HRTF_fname);
    [HRTF] = MTLreadDIR(HRTF);
    
    % THIS CONVERSION IS NOT NECESSARY WITH ANY OF THE BIRDS WE USE NOW:
    % %convert coordinates to double polar, only if not using 929 or 930 or ones
    % if (~strcmp(PDR.HRTF_fname(1:3),'930') & ~strcmp(PDR.HRTF_fname(1:3),'929'))
    %     PDR.HRTF_dir_matrix = sphere2double(PDR.HRTF_dir_matrix);
    % end

    direc = HRTF.dir_matrix; % NOTE: first row = Elevation and 2nd row = Azimuth !!!
    for i=1:PDR.TEST_nlocs
        idx{i}=find(direc(1,:)==PDR.TEST_locs(i,1) & direc(2,:)==PDR.TEST_locs(i,2));
        HRTF.TestL(:,i) = MTLreadCH(idx{i}*2-1, HRTF);
        HRTF.TestR(:,i) = MTLreadCH(idx{i}*2, HRTF);
    end
    
    if PDR.flag_adapt
        % Make FIR coefficients for Adaptor:
        PDR.ADAPT_coefs = makeGammaFIR(PDR.stim_Fs,PDR.ADAPT_cF,PDR.ADAPT_species);
        idxm=find(direc(1,:)==PDR.ADAPT_loc(1,1) & direc(2,:)==PDR.ADAPT_loc(1,2));
        % HRTFs for adapting stimulus:
        HRTF.AdaptL(1,:) = MTLreadCH(idxm*2-1, HRTF);
        HRTF.AdaptR(1,:) = MTLreadCH(idxm*2, HRTF);
    end
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
session.hPlot=figure; whitebg(gcf,'k');
session.hTracePlot=subplot(2,2,1:2);
session.hInfo=subplot(2,2,3); axis off;
session.hTrialPlot=subplot(2,2,4); axis off;
set(session.hPlot,'renderer','OpenGL'); %use OpenGL for renderer
set(session.hPlot, 'Position', [0.05*screen_size(3) 0.05*screen_size(4) 0.7*screen_size(3) 0.8*screen_size(4)] );
txt(1) = text(.01,.9,''); 
txt(2) = text(.01,.7,''); 
txt(3) = text(.01,.5,''); 
txt(4) = text(.01,.3,'');
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

cd(PDR.data_path)
p = pwd;
p = [p '\'];
if strcmp(PDR.data_path,p)
else
    warndlg('Something could be wrong with the path setup!');
    return
end

%write header information to file... saving global variables
PDR = orderfields(PDR); % order fields by ASCII dictionary order
save ([PDR.data_path PDR.filename '.mat'], 'passmein','PDR','HRTF');

% ENGAGE the main spatialPDR script:
spatialPDR;

b=clock;
PDR.stoptime(1:2)=b(4:5);
%write header information to file... saving global variables
cd(PDR.data_path)
save ([PDR.data_path PDR.filename '.mat'], 'passmein','PDR','HRTF','TDT');
cd ..\code