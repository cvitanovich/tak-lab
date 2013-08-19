function C=calib_sound
%% CALIBRATION SCRIPT

%% PARAMETERS:
% note: microphones should also be calibrated for SPL vs. Voltage on the same day

% do I want to use HRTFs? earphones or speakers?
% 1 for yes, 0 for no!
C.HRTFs=1;
C.EPHONES=1;
C.SPEAKERS=0;
% what scales to try?
C.scales_2_try_for_cutoffs=fliplr(10.^([0:.09:4.5 log10(32760)]));
% what atten to use?
C.atten=0;
C.stim_Fs=30000;
C.buf_pts=30000; % 1 sec sounds
C.stim_dur = C.buf_pts/C.stim_Fs;
% mic_type = Knowles or SPL meter
C.MIC_TYPE='Knowles';
C.CALIB_TYPE='Intraural';
% HRTF coefficients file:
if(C.HRTFs==1)
    C.HRTF_fname=''; % what file to use for HRTFs???
    C.HRTF_nlines=255;
    C.HRTF_directory='';
    C.HRTF_Left=[]; % left coefs
    C.HRTF_Right=[]; % right coefs
end
% sounds
C.sounds = {'GTONE','OCTAVE','BBN'};
C.ramp = 5; % 5 ms ramp
C.target_rms=0.1; % desired rms amplitude of sound vector
C.GTONE_cF=6000;
C.GTONE_species = 'owl';
C.OCTAVE_bandwidth=[4000 8000];
C.BBN_bandwidth=[2000 11000];
% filename and data path
c=clock; c=round(c);
yr=num2str(c(1)); yr=yr(end-1:end);
mo=num2str(c(2)); day=num2str(c(3));
C.filename=['Calib' C.CALIB_TYPE '_' yr mo day '_Atten' num2str(C.atten)];

%% CHECK POINT
if strcmp(C.MIC_TYPE,'Knowles')
    h=warndlg({'USING KNOWLES MICS FOR INTRAURAL CALIBRATIONS','CHANNEL 0 = LEFT EAR, CHANNEL 1 = RIGHT EAR ... OKAY???'},'warning');
    uiwait(h);
else
    h=warndlg('Calibrating speaker output with a single mic (SPL Meter?)... OKAY???','warning');
    uiwait(h);
end


%% LOAD HRTF COEFFICIENTS:
if(C.HRTFs==1)
    if(strcmp(C.HRTF_fname((end-3):end),'.mat')) % DOT MAT FORMAT
        LT=zeros(1,C.HRTF_nlines);
        RT=LT;
        EL=0; AZ=0;
        [LT, RT] = readHRTFdotMAT(C.HRTF_directory,C.HRTF_fname,EL,AZ);
        HRTF.left=LT;
        HRTF.right=RT;
    else % using another format (e.g. for 930 or 929)
        %read HRTF coefficients files:
        [HRTF] = MTLreadHDR(C.HRTF_directory,C.HRTF_fname);
        [HRTF] = MTLreadDIR(HRTF);
        
        % THIS CONVERSION IS NOT NECESSARY WITH ANY OF THE BIRDS WE USE NOW:
        % %convert coordinates to double polar, only if not using 929 or 930 or ones
        % if (~strcmp(C.HRTF_fname(1:3),'930') & ~strcmp(C.HRTF_fname(1:3),'929'))
        %     C.HRTF_dir_matrix = sphere2double(C.HRTF_dir_matrix);
        % end
        
        direc = HRTF.dir_matrix; % NOTE: first row = Elevation and 2nd row = Azimuth !!!
        idx=find(direc(1,:)==0 & direc(2,:)==0);
        chan_left = MTLreadCH(idx*2-1, HRTF);
        chan_right = MTLreadCH(idx*2, HRTF);
        HRTF.left=chan_left;
        HRTF.right=chan_right;
    end
else
    % pass sounds through without filtering!
    tmp=zeros(1,255);
    tmp(1)=1;
    HRTF.left=tmp;
    HRTF.right=tmp;
end

%% MAKE BAND LIMITED NOISE
bandwidth=C.OCTAVE_bandwidth;
ok=0; cnt=0; seed=1;
while(~ok)
    C.OCTAVE_seed=seed+cnt;
    tmp= makeTest(C.OCTAVE_seed,C.stim_dur*1000,bandwidth(1),bandwidth(2),C.stim_Fs,0);
    tmp = 0.999.*(tmp ./ (max(abs(tmp)))); % NORMALIZE
    if strcmp(C.MIC_TYPE,'Knowles')
        % Convolve with HRTF (ABL Equalized) for frontal location (El,Az)=(0,0)
        C.OCTAVE_left = filter(HRTF.left,1,tmp);
        C.OCTAVE_right = filter(HRTF.right,1,tmp);
        % NORMALIZE TO +/- 1
        mx=max(abs([C.OCTAVE_left C.OCTAVE_right]));
        C.OCTAVE_left=C.OCTAVE_left./mx;
        C.OCTAVE_right=C.OCTAVE_right./mx;
        rms_left = sqrt(mean(C.OCTAVE_left.^2));
        rms_right = sqrt(mean(C.OCTAVE_right.^2));
        rms_avg = (rms_left+rms_right)/2;
    
        C.OCTAVE_left=(C.target_rms/rms_avg).*C.OCTAVE_left; % sets rms ABL to target rms
        C.OCTAVE_right=(C.target_rms/rms_avg).*C.OCTAVE_right; % sets rms ABL to target rms
        if(max(abs([C.OCTAVE_left C.OCTAVE_right]))<1)
            ok=1; % avoids CLIPPING!
            break
        end
    else
        C.OCTAVE_sound = tmp;
        rms = sqrt(mean(C.OCTAVE_sound.^2));
        C.OCTAVE_sound = (C.target_rms/rms).*C.OCTAVE_sound; % set rms equal to target rms
        if(max(abs(C.OCTAVE_sound))<1)
            ok=1; % avoids CLIPPING!
            break
        end
    end
    cnt=cnt+1;
end
if strcmp(C.MIC_TYPE,'Knowles')
    C.OCTAVE_left = rampMySound(C.OCTAVE_left,C.ramp,C.stim_Fs);
    C.OCTAVE_right = rampMySound(C.OCTAVE_right,C.ramp,C.stim_Fs);
else
    C.OCTAVE_sound = rampMySound(C.OCTAVE_sound,C.ramp,C.stim_Fs);
end

%% MAKE BROADBAND NOISE
bandwidth=C.BBN_bandwidth;
ok=0; cnt=0; seed=1;
while(~ok)
    C.BBN_seed=seed+cnt;
    tmp= makeTest(C.BBN_seed,C.stim_dur*1000,bandwidth(1),bandwidth(2),C.stim_Fs,0);
    tmp = 0.999.*(tmp ./ (max(abs(tmp)))); % NORMALIZE
    if strcmp(C.MIC_TYPE,'Knowles')
        % Convolve with HRTF (ABL Equalized) for frontal location (El,Az)=(0,0)
        C.BBN_left = filter(HRTF.left,1,tmp);
        C.BBN_right = filter(HRTF.right,1,tmp);
        % NORMALIZE TO +/- 1
        mx=max(abs([C.BBN_left C.BBN_right]));
        C.BBN_left= C.BBN_left./mx;
        C.BBN_right= C.BBN_right./mx;
        rms_left= sqrt(mean(C.BBN_left.^2));
        rms_right= sqrt(mean(C.BBN_right.^2));
        rms_avg = (rms_left+rms_right)/2;
        C.BBN_left=(C.target_rms/rms_avg).*C.BBN_left; % sets rms ABL to target rms
        C.BBN_right=(C.target_rms/rms_avg).*C.BBN_right; % sets rms ABL to target rms
        if(max(abs([C.BBN_left C.BBN_right]))<1)
            ok=1; % avoids CLIPPING!
            break
        end
    else
        C.BBN_sound = tmp;
        rms = sqrt(mean(C.BBN_sound.^2));
        C.BBN_sound = (C.target_rms/rms).*C.BBN_sound; % set rms equal to target
        if(max(abs(C.BBN_sound)<1))
            ok=1; % avoids CLIPPING!
            break
        end
    end 
    cnt=cnt+1;
end
if strcmp(C.MIC_TYPE,'Knowles')
    C.BBN_left = rampMySound(C.BBN_left,C.ramp,C.stim_Fs);
    C.BBN_right = rampMySound(C.BBN_right,C.ramp,C.stim_Fs);
else
    C.BBN_sound = rampMySound(C.BBN_sound,C.ramp,C.stim_Fs);
end

%% MAKE GAMMATONE
C.GTONE_coefs = makeGammaFIR(C.stim_Fs,C.GTONE_cF,C.GTONE_species);
ok=0; cnt=0; seed=1;
while(~ok)
    C.GTONE_seed=seed+cnt;
    rand('state',C.GTONE_seed);
    tmp=rand(1,C.buf_pts);
    tmp=filtfilt(C.GTONE_coefs,1,tmp);
    tmp = 0.999.*(tmp ./ (max(abs(tmp)))); % NORMALIZE
    if strcmp(C.MIC_TYPE,'Knowles')
        % Convolve with HRTFs (ABL Equalized) for frontal location (El,Az)=(0,0)
        C.GTONE_left = filter(HRTF.left,1,tmp);
        C.GTONE_right = filter(HRTF.right,1,tmp);
        % NORMALIZE TO +/- 1
        mx=max(abs([C.GTONE_left C.GTONE_right]));
        C.GTONE_left= C.GTONE_left./mx;
        C.GTONE_right= C.GTONE_right./mx;
        rms_left= sqrt(mean(C.GTONE_left.^2));
        rms_right= sqrt(mean(C.GTONE_right.^2));
        rms_avg = (rms_left+rms_right)/2;
        C.GTONE_left=(C.target_rms/rms_avg).*C.GTONE_left; % sets rms ABL to target rms
        C.GTONE_right=(C.target_rms/rms_avg).*C.GTONE_right; % sets rms ABL to target rms
        if(max(abs([C.GTONE_left C.GTONE_right]))<1)
            ok=1; % avoids CLIPPING!
            break
        end
    else
        C.GTONE_sound=tmp;
        rms=sqrt(mean(C.GTONE_sound.^2));
        C.GTONE_sound = (C.target_rms/rms).*C.GTONE_sound; % set rms equal to target
        if(max(abs(C.BBN_sound)<1))
            ok=1; % avoids CLIPPING!
            break
        end
    end
    cnt=cnt+1;
end
if strcmp(C.MIC_TYPE,'Knowles')
    C.GTONE_left = rampMySound(C.GTONE_left,C.ramp,C.stim_Fs);
    C.GTONE_right = rampMySound(C.GTONE_right,C.ramp,C.stim_Fs);
else
    C.GTONE_sound = rampMySound(C.GTONE_sound,C.ramp, C.stim_Fs);
end

keyboard; return;

%% INITIALIZE TDT:
if strcmp(C.MIC_TYPE,'Knowles')
    TDT.nPlayChannels=2;
    TDT.nRecChannels=2; % one for each mic!
    TDT.playpts={C.buf_pts,C.buf_pts};
elseif strcmp(C.MIC_TYPE,'SPL Meter')
    TDT.nPlayChannels=1;
    TDT.nRecChannels=1;
    TDT.playpts={C.buf_pts};
else
    h=warndlg('Improper option for microphone type!','warning');
    uiwait(h);
    return;
end

TDT.recpts=TDT.playpts;
TDT.din = 1;
TDT.Fs = C.stim_Fs;
TDT.npts_total_play=C.buf_pts;
TDT.dec_factor=0; % don't decimate
TDT.srate=1e6 / TDT.Fs;
TDT.buf_pts=C.buf_pts;

% TDT INIT
out=TDT_init;
if(out==-1); return; end;

% PD1 SETUP
if strcmp(C.MIC_TYPE,'Knowles')
    PD1_Route_IntrauralKnowles(TDT);
elseif strcmp(C.MIC_TYPE,'SPL Meter')
    PD1_Route_Speaker(TDT);
else
    h=warndlg('Improper option for microphone type!','warning');
    uiwait(h);
    return;
end

% attenuator settings
TDT_attens(TDT);

%% INITIALIZE BUFFERS

% INITIALIZE BUFFERS
TDT=TDT_buffers(TDT);

if strcmp(C.MIC_TYPE,'Knowles')
    C.TMP_LT=TDT.n_total_buffers+1;
    S232('allotf',C.TMP_LT,TDT.buf_pts);
    C.TMP_RT=TDT.n_total_buffers+2;
    S232('allotf',C.TMP_RT,TDT.buf_pts);
    TDT.n_total_buffers=TDT.n_total_buffers+2;
elseif strcmp(C.MIC_TYPE,'SPL Meter')
    C.TMP=TDT.n_total_buffers+1;
    S232('allotf',C.TMP,TDT.buf_pts);
    TDT.n_total_buffers=TDT.n_total_buffers+1;
end
%% MAIN LOOP:

for cnt=1:length(C.sounds)
    if strcmp(C.MIC_TYPE,'Knowles')
        eval(['left_snd = C.' C.sounds{cnt} '_left;']);
        eval(['right_snd = C.' C.sounds{cnt} '_right;']);
        [lo, hi] = det_cutoffs(C.scales_2_try_for_cutoffs, left_snd, right_snd);
    else % using SPL Meter (speaker)
        eval(['snd = C.' C.sounds{cnt} '_sound;']);
        [lo, hi] = det_cutoffs(C.scales_2_try_for_cutoffs, snd);
    end
	lo=round(10^lo); hi=round(10^hi);
    scales=lo:((hi-lo)/C.nscales):hi;
    if strcmp(C.MIC_TYPE,'Knowles')
        rms = test_scales(scales, left_snd, right_snd);
    else
        rms = test_scales(scales, snd);
    end
    dbs = (1/C.coeffs(2))*log( (rms - C.coeffs(3)) ./ C.coeffs(1) );
    
    % REGRESSION FIT
    h=figure;
    hold on;
    xrange = log10(1:1:32760);
    xes = log10(scales);
    yes = dbs;
    colr = [1 0 0];
    [rsq, coefs] = regress_stats(xes,yes,0.05,xrange,colr,1);
    eval(['C.' C.sounds{cnt} '_RSQ = rsq;']);
    eval(['C.' C.sounds{cnt} '_COEFS= coefs;']);
    % title string
    title_string=[];
    title_string{1} = ['SOUND IS: ' C.sounds{cnt}];
    title_string{2} = ['RSQ = ' num2str(rsq)];
    title(title_string,'FontSize',8);
    % axes labels
    xlabel('log10(scales)','FontSize',8);
    if strcmp(C.MIC_TYPE,'Knowles')
        ylabel('SPL (dB, ABL)','FontSize',8);
    else
        ylabel('SPL (dB)','FontSize',8);
    end
    hold off;
    % save this figure
    set(h,'InvertHardcopy','off');
    fname = [C.data_path C.filename C.sounds{cnt} '_fit'];
    saveas(h,fname,'fig');
end
%write header information to file... saving global variables
save ([C.data_path C.filename '.mat'],'C','TDT','HRTF');
str{1} = 'Variables saved in: ';
str{2} = [C.data_path C.filename '.mat'];
hMsg=msgbox(str); uiwait(hMsg);


%% SUBROUTINES
function [lo, hi] = det_cutoffs(scales2try,snd0,snd1)
global MIC
if(nargin==3)
    left_snd=snd0; clear(snd0);
    right_snd=snd1; clear(snd1);
    rms = test_scales(scales2try,left_snd,right_snd);
else
    snd=snd0; clear(snd0);
    rms = test_scales(scales2try,snd);
end
% select cutoffs
dbs = (1/MIC.coeffs(2))*log( (rms - MIC.coeffs(3)) ./ MIC.coeffs(1) );
[lo, hi] = select_cutoffs(log10(scales2try),dbs);
if hi>log10(32760)
    hi=log10(32760);
end
if lo<0
    lo=0;
end

function [lo, hi] = select_cutoffs(xes,yes)
% Plots data in a scatter plot and asks user to input lo/hi cuttoffs
% for fitting linear portion of data range
screen_size = get(0, 'ScreenSize');
hFig=figure;
set(hFig, 'Position', [0.02*screen_size(3) 0.05*screen_size(4) 0.8*screen_size(3) 0.8*screen_size(4)] );

scatter(xes,yes,'.');
prompt = {'Enter Lower Cutoff','Enter Upper Cutoff'};
dlg_title = 'Input Cutoffs:';
num_lines = 1;
def = {'',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
lo=str2num(answer{1});
hi=str2num(answer{2});
close(hFig);
drawnow;

function rms = test_scales(scales,snd0,snd1)
global TDT C
if(nargin==3)
    % two outputs (headphones)
    left_snd=snd0; clear(snd0);
    right_snd=snd1; clear(snd1);
    s232('pushf',left_snd, TDT.buf_pts);
    s232('qpopf',C.TMP_LT);
    s232('pushf',right_snd, TDT.buf_pts);
    s232('qpopf',C.TMP_RT);
else
    % one output (speaker)
    snd=snd0; clear(snd0);
    s232('pushf',snd,TDT.buf_pts);
    s232('qpopf',C.TMP)
end
 hWait=waitbar(0,'playing sounds...');

for j=1:length(scales)
    if(nargin==3)
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
    if(nargin==3)
        S232('qpush16', TDT.rec_buffers{1}(1));
        left_tmp=S232('pop16');
        S232('qpush16', TDT.rec_buffers{2}(1));
        right_tmp=S232('pop16');
        % calculate rms voltage
        left_rms = sqrt(mean(left_tmp.^2));
        right_rms = sqrt(mean(right_tmp.^2));
        rms(j) = (left_rms + right_rms)/2;
    else
        s232('qpush16',TDT.rec_buffers{1}(1));
        rms=s232('pop16'); % one mic (no need to average)
    end 
	waitbar(j/length(scales),hWait);
end
close(hWait);

% REGRESSION STATS
function [RSQ, COEFFS] = regress_stats(xes,yes,alfa,xrange,col,info_flg)
% plot data
n=length(yes);
df=length(yes)-2;
for j = 1:n
    plot(xes,yes,'.','Color',col);
end
% determine fit coefficients and CI
[p, S] = polyfit(xes,yes,1);
[Y,DELTA]=polyconf(p,xrange,S,alfa);
m=p(1);
b=p(2);
yfit =  m.* xrange + b;
ypred = (m.*xes + b);
yresid = yes - ypred;
SSresid = sum(yresid.^2);
SStotal = (length(yes)-1) * var(yes);
RSQ = 1 - SSresid/SStotal;
% plot fit with CI
hold on;
plot(xrange,Y,'c-','Color',col);
plot(xrange,Y-DELTA,'r--','Color',col);
plot(xrange,Y+DELTA,'r--','Color',col);
hold off;
% standard error calculation
SSxx = sum(xes.^2)-n*mean(xes)^2;
SSyy = sum(yes.^2)-n*mean(yes)^2;
SSxy = sum(xes.*yes)-n*mean(xes)*mean(yes);
s=sqrt((SSyy-m*SSxy)/(n-2));
SE_m = s/sqrt(SSxx); % standard error for slope
SE_b = s*sqrt((1/n)+mean(xes)^2/SSxx); % standard error for intercept
% determine t statistic
step = 0.01;
t=step;
cum=0;
while cum < (1-alfa)
    tes=-t:step:t;
    tmp=tpdf(tes,df);
    cum=sum(tmp)*step;
    if cum > 0.95
        break;
    else
        t=t+step;
    end
end
% determine coefficient CIs
COEFFS = zeros(3,3);
COEFFS = [m m-SE_m*t m+SE_m*t; b b-SE_b*t b+SE_b*t];
% set axis values
minX=min(xrange);
maxX=max(xrange);
minY=min(Y-DELTA);
maxY=max(Y+DELTA);
axis([minX maxX minY maxY]);
if info_flg
    % plot info
    eqn = ['Fit: y = mx+b = ' num2str(m) 'x + ' num2str(b)];
    text(minX+0.1*(maxX-minX),minY+0.9*(maxY-minY),eqn,'FontSize',8);
    mcoeff=['m = ' num2str(m) ' [' num2str(COEFFS(1,2)) ',' num2str(COEFFS(1,3)) ']'];
    text(minX+0.1*(maxX-minX),minY+0.75*(maxY-minY),mcoeff,'FontSize',8);
    bcoeff=['b = ' num2str(b) ' [' num2str(COEFFS(2,2)) ',' num2str(COEFFS(2,3)) ']'];
    text(minX+0.1*(maxX-minX),minY+0.6*(maxY-minY),bcoeff,'FontSize',8);
end

function PD1_Route_IntrauralKnowles(TDT)
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
s232('PD1specOB',TDT.din,s232('OB',0),s232('ADC',0));
% LED thresholds
S232('PD1setIO',TDT.din,0.01,9.99,0.01,9.99);

function load_HRTF_coefs(HRTF)
global C TDT
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
s232('pushf',HRTF.TestL,C.HRTF_nlines);
s232('PreLoadRaw',TDT.din,s232('DSPid',0),MONO,STACK,[],[],1,1,1);
s232('pushf',HRTF.TestR,C.HRTF_nlines);
s232('PreLoadRaw',TDT.din,s232('DSPid',1),MONO,STACK,[],[],1,1,1);
    