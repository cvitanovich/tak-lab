function C=calib_sound(C,HRTF)
%% CHECK POINT
if strcmp(C.MIC_TYPE,'Knowles')
    h=warndlg({'USING KNOWLES MICS FOR INTRAURAL CALIBRATIONS','CHANNEL 0 = LEFT EAR, CHANNEL 1 = RIGHT EAR ... OKAY???'},'warning');
    uiwait(h);
else
    h=warndlg('Calibrating speaker output with a single mic (SPL Meter?)... OKAY???','warning');
    uiwait(h);
end
if(C.EPHONES)
    h=warndlg('Using earphones, not speaker, okay?','warning');
    uiwait(h);
else
    h=warndlg('Using speaker, okay?','warning');
    uiwait(h);
end

%% INITIALIZE TDT:
if(C.EPHONES)
    TDT.nPlayChannels=2;
    TDT.playpts={C.buf_pts,C.buf_pts};
else
    TDT.nPlayChannels=1;
    TDT.playpts={C.buf_pts};
end

if strcmp(C.MIC_TYPE,'Knowles')
    TDT.nRecChannels=2; % one for each mic!
    TDT.recpts={C.buf_pts,C.buf_pts};
elseif strcmp(C.MIC_TYPE,'SPL Meter')
    TDT.nRecChannels=1;
    TDT.recpts={C.buf_pts};
else
    h=warndlg('Improper option for microphone type!','warning');
    uiwait(h);
    return;
end

TDT.din = 1;
TDT.Fs = C.stim_Fs;
TDT.npts_total_play=C.buf_pts;
TDT.dec_factor=0; % don't decimate
TDT.srate=1e6 / TDT.Fs;
TDT.buf_pts=C.buf_pts;

% TDT INIT
out=TDT_init;
if(out==-1)
    h=warndlg('RESET S232 Controller, okay?','warning'); uiwait(h);
end

% PD1 SETUP
if(C.EPHONES)
    PD1_Route_Earphones(TDT);
else
    PD1_Route_Speaker(TDT);
end

%% INITIALIZE BUFFERS
TDT=TDT_buffers(TDT);
if(C.EPHONES)
    C.TMP_LT=TDT.n_total_buffers+1;
    S232('allotf',C.TMP_LT,TDT.buf_pts);
    C.TMP_RT=TDT.n_total_buffers+2;
    S232('allotf',C.TMP_RT,TDT.buf_pts);
    TDT.n_total_buffers=TDT.n_total_buffers+2;
else
    C.TMP=TDT.n_total_buffers+1;
    S232('allotf',C.TMP,TDT.buf_pts);
    TDT.n_total_buffers=TDT.n_total_buffers+1;
end

%% LOAD HRTFs IF USED
if(C.EPHONES && C.HRTFs)
    load_HRTF_coefs(HRTF,TDT);
elseif(C.EPHONES && ~C.HRTFs)
    % ephones w/o filtering
    HRTF.nlines = 255;
    HRTF.left = [0.99999 zeros(1,254)];
    HRTF.right = HRTF.left;
else
    % do nothing
end

% date
c=clock; c=round(c);
yr=num2str(c(1)); yr=yr(end-1:end);
mo=num2str(c(2)); day=num2str(c(3));

for j=1:length(C.sounds)
    % sound type
    current_snd=C.sounds{j};
    
    % filename
    C.filename=['Calib' C.CALIB_TYPE C.MIC_TYPE];
    if(C.HRTFs)
        C.filename=[C.filename '_HRTFs'];
    end
    if(C.EPHONES)
        C.filename=[C.filename '_EPHONES'];
    end
    if(C.SPEAKERS)
        C.filename=[C.filename '_SPEAKERS'];
    end
    C.filename=[C.filename '_' current_snd '_' yr mo day];
    
    ATTENS=[];SCALES=[];VOLTAGES=[];
    for k=1:length(C.attens)
        current_atten=C.attens(k);
        % attenuator settings
        for(j=1:TDT.nPlayChannels)
            S232('PA4atten',j,current_atten);
        end
        ATTENS=[ATTENS current_atten.*ones(1,length(C.scales))];
        SCALES=[SCALES C.scales];
        %% PLAY SOUNDS AND RECORD VOLTAGE FROM MICS
        if(C.EPHONES)
            eval(['left_snd = C.' current_snd '_sound;']);
            eval(['right_snd = C.' current_snd '_sound;']);
            v = test_scales(C,TDT,C.scales, left_snd, right_snd);
        else % using SPL Meter (speaker)
            eval(['snd = C.' current_snd '_sound;']);
            v = test_scales(C,TDT,C.scales, snd);
        end
        VOLTAGES=[VOLTAGES v];
        close all; figure; plot(C.scales,v);
        xlabel('SCALES'); ylabel('RMS VOLTAGE (mV)'); title([current_snd ', Atten = ' num2str(current_atten)]);
    end
    
    %write data to a text file:
    fname=[C.data_path C.filename 'A'];
    count = double('A'+0);
    while exist ([fname '.txt'],'file');
    count = count + 1;
        if count > 90
            disp('There are already several files with similar names!');
            fname = input('Enter a unique filename for this data file: ', 's');
            break;
        else
        fname(end) = char(count);
        end
    end
    fid = fopen([fname '.txt'],'w');
    DATA=[ATTENS; SCALES; VOLTAGES];
    formatSpec='%3.1f\t%5.0f\t%5.3f\n';
    nbytes=fprintf(fid,formatSpec,DATA);
    fclose(fid);
    disp(['Saved data as: ' C.filename '.txt']);
    disp(['Wrote ' num2str(nbytes/1000) ' kBytes of calibration data...']);
    
    %write header information to file... saving global variables
    if(C.HRTFs)
        save ([fname '.mat'],'C','TDT','HRTF');
    else
        save ([fname '.mat'],'C','TDT');
    end
    disp(['Variables saved in: ' C.data_path C.filename '.mat']);
end

%% SUBROUTINES
function rms = test_scales(C,TDT,scales,current_snd0,current_snd1)
% filename for saving voltage trace
FNAME=[C.data_path C.filename];
left_rms=[]; right_rms=[];
if(nargin>4)
    % two outputs (headphones)
    left_current_snd=current_snd0;
    right_current_snd=current_snd1;
    s232('pushf',left_current_snd, TDT.buf_pts);
    s232('qpopf',C.TMP_LT);
    s232('pushf',right_current_snd, TDT.buf_pts);
    s232('qpopf',C.TMP_RT);
else
    % one output (speaker)
    current_snd=current_snd0;
    s232('pushf',current_snd,TDT.buf_pts);
    s232('qpopf',C.TMP)
end
hWait=waitbar(0,'playing sounds...');
for j=1:length(scales)
    s232('dropall');
    if(nargin>4)
        % two outputs (headphones)
        S232('qpushf',C.TMP_LT);
        S232('scale',scales(j));
        s232('qpop16',TDT.stim_buffers{1}(1));
        s232('qpushf',C.TMP_RT);
        S232('scale',scales(j));
        s232('qpop16',TDT.stim_buffers{2}(1));
    else
        % one output (speaker)
        S232('qpushf',C.TMP);
        S232('scale',scales(j));
        s232('qpop16',TDT.stim_buffers{1}(1));
    end
    S232('seqplay',TDT.play_spec);
    % recording voltage
    s232('seqrecord',TDT.rec_spec);
    s232('PD1arm',1);
    s232('PD1go',1);
    while(s232('PD1status',1))
        % wait for TDT to finish
    end
    % pop recorded voltage traces into matlab
    if exist('left_tmp'); clear left_tmp; end;
    if exist('right_tmp'); clear right_tmp; end;
    if exist('tmp'); clear tmp; end;
    if(strcmp(C.MIC_TYPE,'Knowles'))
        % get record buffer for left ear
        S232('qpush16', TDT.rec_buffers{1}(1));
        S232('qdup');
        left_tmp=S232('pop16');
        % save to disc
        S232('make',0,1);
        S232('qpop16',TDT.rec_buffers{1}(1));
        S232('dama2disk16',TDT.rec_buffers{1}(1),[FNAME '_left.rec'],1);
        % get record buffer for right ear
        S232('qpush16', TDT.rec_buffers{2}(1));
        S232('qdup');
        right_tmp=S232('pop16');
        % save to disc
        S232('make',0,1);
        S232('qpop16',TDT.rec_buffers{2}(1));
        S232('dama2disk16',TDT.rec_buffers{2}(1),[FNAME '_right.rec'],1);
        % calculate rms voltage
        left_tmp = left_tmp - mean(left_tmp);
        right_tmp = right_tmp - mean(right_tmp);
        left_rms(end+1) = sqrt(mean(left_tmp.^2));
        right_rms(end+1) = sqrt(mean(right_tmp.^2));
        rms(j) = (left_rms(end) + right_rms(end))/2;
    else
        s232('qpush16',TDT.rec_buffers{1}(1));
        S232('qdup');
        tmp=s232('pop16'); % one mic (no need to average)
        % save to disc
        S232('make',0,1);
        S232('qpop16',TDT.rec_buffers{1}(1));
        S232('dama2disk16',TDT.rec_buffers{1}(1),[FNAME '.rec'],1);
        tmp = tmp - mean(tmp);
        rms(j)=sqrt(mean(tmp.^2));
    end 
	waitbar(j/length(scales),hWait);
end
close(hWait);

function PD1_Route_Earphones(TDT)
% routing for knowles recording (intraural)
% commands to set up routing in the TDT
SRATE = 1e6 / TDT.Fs;
S232('PD1clear',TDT.din);
S232('PD1fixbug',TDT.din);
S232('PD1nstrms',TDT.din,2,2);
S232('PD1srate',TDT.din,SRATE);
S232('PD1npts',TDT.din,TDT.npts_total_play);
S232('PD1clrsched',TDT.din);
% left earphone routing
s232('PD1addsimp',TDT.din,s232('DSPout',0),s232('DAC',0));
% right earphone routing
s232('PD1addsimp',TDT.din,s232('DSPout',1),s232('DAC',1));
% IREGS
S232('PD1addsimp',TDT.din,s232('IREG',0),s232('DSPin',0));
S232('PD1addsimp',TDT.din,s232('IREG',0),s232('DSPin',1));
% sound routing
s232('PD1specIB',TDT.din,s232('IB',0),s232('IREG',0));
% recording (two channels: left/right knowles!)
s232('PD1specOB',TDT.din,s232('OB',0),s232('ADC',0));
s232('PD1specOB',TDT.din,s232('OB',1),s232('ADC',1));
% LED thresholds
S232('PD1setIO',TDT.din,0.01,9.99,0.01,9.99);

function PD1_Route_Speaker(TDT)
% routing for just a simple 1 speaker output!
% commands to set up routing in the TDT
SRATE = 1e6 / TDT.Fs;
S232('PD1clear',TDT.din);
S232('PD1fixbug',TDT.din);
S232('PD1nstrms',TDT.din,2,2);
S232('PD1srate',TDT.din,SRATE);
S232('PD1npts',TDT.din,TDT.npts_total_play);
S232('PD1clrsched',TDT.din);
% sound routing
s232('PD1specIB',TDT.din,s232('IB',0),s232('DAC',0));
% recording
for(j=0:(TDT.nRecChannels-1))
    s232('PD1specOB',TDT.din,s232('OB',j),s232('ADC',j));
end
% LED thresholds
S232('PD1setIO',TDT.din,0.01,9.99,0.01,9.99);

function load_HRTF_coefs(HRTF,TDT)
% opmode parameters
MONO=1;
STEREO=2;
MONSTER=3;
% stype params
DAMA_F=5;
STACK=6;
% DSP IDs:
DSP_Test=0;
DSP_Adapt=1;
% clear stack
s232('dropall');
% idle ALL dsps
S232('PD1idleDSP',TDT.din,hex2dec('0FFFFFFF'));
% load left/right coefs to DSP
s232('pushf',HRTF.left,HRTF.nlines);
s232('PreLoadRaw',TDT.din,s232('DSPid',0),MONO,STACK,[],[],1,1,1);
s232('pushf',HRTF.right,HRTF.nlines);
s232('PreLoadRaw',TDT.din,s232('DSPid',1),MONO,STACK,[],[],1,1,1);
    