function KNOWLES=calib_knowles
% changed to require a single speaker... 9/21/2013
% STIMULUS SETUP
KNOWLES = struct(...
    'minfreq',2000,...
    'maxfreq',11000,...
    'Fs',30000,...
    'duration',5,...
    'ramp',5,...
    'rms',0.1,...
    'attens',30:-5:10,...
    'scales',3000:3000:30000,...
    'data_path','',...
    'filename',[]);

if(~ispc && strcmp(getenv('USER'),'cvitanovich'))
    data_path='/Users/cvitanovich/Documents/MATLAB/data/';
else
    data_path='E:\andrew\data\calib\';
end

KNOWLES.data_path=data_path;

% setup filename
c=clock;
KNOWLES.filename=['knowles_calib_' num2str(c(1)) '_' num2str(c(2)) '_' num2str(c(3)) '_A'];
count = double('A'+0);
while exist ([KNOWLES.data_path KNOWLES.filename '.mat'],'file');
    count = count + 1;
    if count > 90
        disp('There are already several files with similar names!');
        KNOWLES.filename = input('Enter a unique filename for this session: ', 's');
        break;
    else
        KNOWLES.filename(end) = char(count);
    end
end

[X,Xtime] = whnoise(KNOWLES.minfreq,KNOWLES.maxfreq,KNOWLES.Fs,KNOWLES.duration);
KNOWLES.buf_pts = length(Xtime);
test=Xtime;
rms = sqrt(mean(test.^2));
test=(KNOWLES.rms/rms).*test;
% ramp the test sound (eliminates speaker popping)
% SOUNDS_ramp should be in ms
ramp_pts = ceil(KNOWLES.Fs*(KNOWLES.ramp/1000));
npts = length(test) - 2*ramp_pts - 2;
rampenv = [0:(1/ramp_pts):1 ones(1,npts) 1:-(1/ramp_pts):0];
test = test' .* rampenv;
% store test sound:
KNOWLES.test = test;
% INITIALIZE TDT
% initialize all SystemII hardware
if(S232('S2init',0,'INIT_PRIMARY',1000)==0)
    warndlg('FAILED to initialize a primary process','warning');
    return;
else
    disp('Initialized primary process');
end

if(S232('APlock',100, 0)==0)
    warndlg('FAILED to get AP lock for initialization','warning');
    return;
else
    disp('AP lock in initiation');
end
if(S232('XBlock',100, 0)==0)
    warndlg('FAILED to get XB lock for initialization','warning');
    S232('APunlock',0);
    warndlg('AP unlock in initiation','warning');
    return;
else
    disp('XB lock in initiation');
end

err = S232('getS2err');
switch err
    case 0
    case 1: warndlg('APOS error in initiation','warning'); return;
    case 2: warndlg('XBUS error in initiation','warning'); return;
end


SRATE =  (1/KNOWLES.Fs)*10^6; % sampling rate for TDT

% PREPARE PD1 FOR CONVERSION:
% one speaker output, two knowles inputs
s232('PD1clear',1);
s232('PD1fixbug',1);
s232('PD1nstrms',1, 1, 2);
s232('PD1srate',1,SRATE);
s232('PD1npts',1, KNOWLES.buf_pts);
s232('PD1clrsched',1);
s232('PD1addsimp',1, s232('IREG',0), s232('DAC',0));
s232('PD1specIB',1, s232('IB',0), s232('IREG',0));
s232('PD1setIO',1,.01,9.99,.01,9.99);
s232('PD1specOB',1,s232('OB',0), s232('ADC',0));
s232('PD1specOB',1,s232('OB',1), s232('ADC',1));
% SET PA4 ATTENUATION:
s232('PA4atten',1,0); % NO ATTEN YET (SET LATER)
s232('PA4atten',2,0); % NO ATTEN YET (SET LATER)

% allot buffers for voltage trace record
if s232('APlock',100,0)==1
    NPTS = KNOWLES.buf_pts;
    
    % play spec list
    CHA_SEQ = s232('_allot16', 10);
    PLAY_SPEC = s232('_allot16',10);
    s232('dpush',10);
    s232('value',0);
    s232('make',0,CHA_SEQ);
    s232('make',1,0);
    s232('qpop16',PLAY_SPEC);
    
    PLAYBUF = s232('_allot16', NPTS);
    
    BUF_TEST = s232('_allotf', NPTS);
    s232('pushf',KNOWLES.test, NPTS);
    s232('qpopf',BUF_TEST);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % play sequences
    s232('dpush', 10);
    s232('value', 0);
    s232('make', 0, PLAYBUF);
    s232('make', 1, 1);
    s232('make', 2, 0);
    s232('qpop16', CHA_SEQ);
    
    % for recording Knowles mic voltage trace
    RECBUF_L1 = s232('_allot16', NPTS);
    RECBUF_R1 = s232('_allot16', NPTS);
    
    % record spec list
    REC_SPEC = s232('_allot16', 10);
    REC_LEFT_SEQ = s232('_allot16', 10);
    REC_RIGHT_SEQ = s232('_allot16', 10);
    
    s232('dpush',10);
    s232('value',0);
    s232('make',0,REC_LEFT_SEQ);
    s232('make',1,REC_RIGHT_SEQ);
    s232('make',2,0);
    s232('qpop16', REC_SPEC);
    
    % record sequences
    s232('dpush',10);
    s232('value',0);
    s232('make',0, RECBUF_L1);
    s232('make',1,1);
    s232('make',2,0);
    s232('qpop16', REC_LEFT_SEQ);
    
    s232('dpush',10);
    s232('value',0);
    s232('make',0, RECBUF_R1);
    s232('make',1,1);
    s232('make',2,0);
    s232('qpop16', REC_RIGHT_SEQ);
    
    s232('APunlock', 0);
else
    disp('Cannot allocate buffers... quitting!')
    return;
end

scalez=KNOWLES.scales;
attens=KNOWLES.attens;
% loop through a range of scales and get SPL readings

step=10; cntdwn=20; cnt=0; spls=[];
str{1}=['Get ready to measure SPLs for ' num2str(length(scalez)*length(attens)) ' sounds'];
str{2}=['in steps of ' num2str(step)];
hMsg=msgbox(str); 
uiwait(hMsg);
% countdown
hCnt=figure; axis off;
for secs=cntdwn:-1:0
    hTxt=text(.5,.5,['Next test in: ' num2str(secs) ' seconds!'],'HorizontalAlignment','center','FontSize',30);
    pause(1);
    delete(hTxt);
end
close(hCnt)
hWait = waitbar(0,'Testing Knowles...');
scalez=scalez(randperm(length(scalez)));
voltages=NaN*ones(1,length(scalez));

% initialize variables
KNOWLES.attens=[];
KNOWLES.scales=[];
KNOWLES.voltages=[];
KNOWLES.spls=[];
num_trials=length(attens)*length(scalez);
prog=0; left_rms=[]; right_rms=[];
for hh=1:length(attens)
    cnt=0; spls=[]; voltages=[];
    % set attens
    s232('PA4atten',1,attens(hh));
    s232('PA4atten',2,attens(hh));
    ii=0;
    while ii<length(scalez)
        ii=ii+1;
        cnt=cnt+1;
        prog=prog+1;
        scale=scalez(ii);
        pause(1); % pause 
        s232('dropall'); % clear the stack
        % play sound
        s232('qpushf',BUF_TEST);
        s232('scale',scale);
        s232('qpop16',PLAYBUF);
        s232('seqplay',PLAY_SPEC);
        % recording voltage
        s232('seqrecord',REC_SPEC);
        s232('PD1arm',1);
        s232('PD1go',1);
        while(s232('PD1status',1))
            % wait for TDT to finish
        end
        % pop recorded voltage traces into matlab
        if exist('left_tmp')
            clear left_tmp;
        end
        if exist('right_tmp')
            clear right_tmp;
        end
        s232('qpush16', RECBUF_L1);
        left_tmp=s232('pop16');
        s232('qpush16', RECBUF_R1);
        right_tmp=s232('pop16');
        % calculate rms voltages and average across both ears
        left_tmp = left_tmp - mean(left_tmp);
        right_tmp = right_tmp - mean(right_tmp);
        % calculate rms voltages
        left_rms(end+1) = sqrt(sum(left_tmp.^2)./length(left_tmp));
        right_rms(end+1) = sqrt(sum(right_tmp.^2)./length(right_tmp));
        avg_rms = (left_rms(end) + right_rms(end))/2;
        voltages(ii)=avg_rms;
        def{cnt}='';
        prompt{cnt}=['Trial #' num2str(prog)];
        % get spl readings if step reached

        if cnt == step || ii==length(scalez)
            button = questdlg('Do you want to repeat those trials?',...
                'Continue Operation','Yes','No','No');
            if strcmp(button,'No')
                if ii==length(scalez) && rem(length(scalez),step)
                    inc = rem(length(scalez),step);
                else
                    inc = step;
                end
                answer=inputdlg(prompt,['Enter SPL Readings for trials ' num2str(prog-inc+1) ' through ' num2str(prog)],1,def);
                for i0=1:size(answer,1)
                    spls(ii-inc+i0)=str2num(answer{i0});
                end
            end
            
            % pause for user & provide countdown
            if prog~=num_trials
                if strcmp(button,'No')
                    numbr=step;
                    if rem(length(scalez),step) && ii~=length(scalez)
                        numbr=rem(length(scalez),step);
                    end
                    str=['Ready for next ' num2str(numbr) ' sounds???'];
                else
                    str=['Repeating last ' num2str(cnt) ' trials again...'];
                    prog = prog - cnt;
                    ii=ii-cnt;
                end
                hMsg=msgbox(str); 
                uiwait(hMsg);
                hCnt=figure; axis off;
                for secs=cntdwn:-1:0
                    hTxt=text(.5,.5,['Next test in: ' num2str(secs) ' seconds!'],'HorizontalAlignment','center','FontSize',30);
                    pause(1)
                    delete(hTxt);
                end
                close(hCnt)
            end
            clear def prompt
            cnt=0;
        end
        waitbar((prog/num_trials),hWait);
    end
    % store sorted spls/scales in cell array
    scalez_sorted=sort(scalez);
    for id=1:length(scalez)
        spls_sorted(id)=spls(find(scalez==scalez_sorted(id)));
        voltages_sorted(id)=voltages(find(scalez==scalez_sorted(id)));
    end
    KNOWLES.attens=[KNOWLES.attens attens(hh)*ones(1,length(scalez))];
    KNOWLES.scales=[KNOWLES.scales scalez_sorted];
    KNOWLES.voltages=[KNOWLES.voltages voltages_sorted];
    KNOWLES.spls=[KNOWLES.spls spls_sorted];
end
close(hWait)
% REGRESSION FIT FOR KNOWLES CALIBRATION
screen_size = get(0, 'ScreenSize');
hKnowles=figure;
set(hKnowles, 'Position', [0.02*screen_size(3) 0.05*screen_size(4) 0.95*screen_size(3) 0.8*screen_size(4)] );
a=0.05; %alpha
hold on;
% exponential regression of avg rms
colr = [0 0 1];
% combine data for all attenuations
xes=KNOWLES.spls;
yes=KNOWLES.voltages;
xind = (max(xes)-min(xes))/100;
xrange = min(xes):xind:max(xes);
beta0=[.05 .2 min(yes)]; % ballpark guess of coefficients for exponential fit
[KNOWLES.rsquared, KNOWLES.coefs] = regress_stats_exponential(xes,yes,a,xrange,colr,beta0);
ax_lead = gca; set(ax_lead,'XColor',colr,'YColor',colr);
% title string
title_string{1} = ['Binaural Mean RMS Voltage vs. Decibels (SPL)'];
title_string{2} = ['RSQ = ' num2str(KNOWLES.rsquared)];
title_string{3} = ['RMS Voltage = ' num2str(KNOWLES.coefs(1)) '*' 'e^{SPL*' num2str(KNOWLES.coefs(2)) '} + ' num2str(KNOWLES.coefs(3))];
title(title_string,'FontSize',8);
% axes labels
xlabel('Level (dB SPL)','FontSize',8);
ylabel('RMS Voltage (averaged across mics)','FontSize',8);
axis tight
hold off;
set(gcf,'InvertHardcopy','off');

fname = [KNOWLES.data_path KNOWLES.filename '_regressfit'];
saveas(hKnowles,fname,'fig');

%write header information to file... saving global variables
save ([KNOWLES.data_path KNOWLES.filename '.mat'],'KNOWLES');
clear str
str{1} = 'Variables saved in: ';
str{2} = [KNOWLES.data_path KNOWLES.filename '.mat'];
hMsg=msgbox(str); uiwait(hMsg);

