function C=calib_sound(C)
%% CHECK POINT
if strcmp(C.MIC_TYPE,'Knowles')
    h=warndlg({'USING KNOWLES MICS FOR INTRAURAL CALIBRATIONS','CHANNEL 0 = LEFT EAR, CHANNEL 1 = RIGHT EAR ... OKAY???'},'warning');
    uiwait(h);
else
    h=warndlg('Calibrating speaker output with a single mic (SPL Meter?)... OKAY???','warning');
    uiwait(h);
end

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
for(j=TDT.nPlayChannels)
    S232('PA4atten',j,C.current_atten);
end

%% INITIALIZE BUFFERS
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

%% PLAY SOUNDS AND RECORD VOLTAGE FROM MICS, THEN ANALYZE
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
fname = [C.data_path C.filename '_fit'];
saveas(h,fname,'fig');
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
    