function spatialPDR
% A FUNCTION TO RUN A SPATIAL PDR EXPERIMENT WITH AN ADAPTOR
global PDR TDT HRTF session

%% INITIATE TDT PARAMETERS HERE
TDT.nPlayChannels=2;
TDT.nBuffers=2;
TDT.bufpts=[PDR.buf_pts PDR.buf_pts]; % buffer points for each play segment
if(PDR.record)
    TDT.nRecChannels=2;
else
    TDT.nRecChannels=0;
end
TDT.dec_factor=PDR.decimationfactor;
TDT.din = 1;
TDT.Fs = PDR.stim_Fs;
TDT.npts_total_play=PDR.npts_totalplay;
TDT.outFN1=[PDR.filename '_REC1.vrt'];
TDT.outFN2=[PDR.filename '_REC2.vrt'];
TDT.ntrials=PDR.ntrials;
TDT.srate=1e6 / TDT.Fs;
TDT.display_flag=1; % flag to display trace during session
TDT.max_signal_jitter=PDR.ISI_jitter;
TDT.disp_proc_time=1; % flag to display processing time for each buffer segment

%% CHECK PARAMETERS
out=check_params;
if(out==-1)
    return;
end

%% INITIALIZE DOUBLE BUFFER PROGRESS VARIABLES
signalScale=0;
readflag=0;
Signalcnt=1; % should be 1 (initialized to 0 in the original C code)
signalScale=0;
record=(TDT.nRecChannels>0);
cnt=1; % ISI counter

%% INITIALIZE TDT
out=TDT_init;
if(out==-1); return; end;

%% INITIALIZE BUFFERS
TDT_buffers;

% stimulus buffer ids
BUF_L1=TDT.stim_buffers(1,1);
BUF_L2=TDT.stim_buffers(1,2);
BUF_R1=TDT.stim_buffers(2,1);
BUF_R2=TDT.stim_buffers(2,2);


%% ADAPTOR FILTERING PARAMS
% initial state
adapt_state = PDR.ADAPT_seed;
% circular buffer for continuous filtered adaptor:
CIRC_BUFS.adaptor=zeros(1,(length(PDR.adapt_coefs)+PDR.buf_pts));
% circular buffers for HRTF filtering (left/right):
CIRC_BUFS.left=circ_buffer; CIRC_BUFS.right=circ_buffer;

%% INITIALIZE PD1
PD1_init;

%% FILTER TEST SOUNDS WITH HRTFS
% buffer assignments for TDT
TESTLEFT=(TDT.n_total_buffers+1):(TDT.n_total_buffers+PDR.TEST_nlocs);
TDT.n_total_buffers=TESTLEFT(end);
TESTRIGHT=(TDT.n_total_buffers+1):(TDT.n_total_buffers+PDR.TEST_nlocs);
TDT.n_total_buffers=TESTRIGHT(end);
filtered_test=zeros(1,PDR.buf_pts);
% loop through test sounds
for(j=1:PDR.TEST_nlocs)
    S232('allotf',PDR.buf_pts);
    S232('allotf',PDR.buf_pts);
    % filter each test sound with HRTFS & Store on AP2 Card
    filtered_test=filter(HRTF.TestL(:,j),1,PDR.TEST_sound);
    S232('pushf',filtered_test,PDR.buf_pts);
    S232('qpopf',TESTLEFT(j));
    filtered_test=filter(HRTF.TestR(:,j),1,PDR.TEST_sound);
    S232('pushf',filtered_test,PDR.buf_pts);
    S232('qpopf',TESTRIGHT(j));
end

%% SET ATTENS
TDT_attens([PDR.base_atten PDR.base_atten]);

%% SIGNAL RAMP BUFFER
S232('allotf',TDT.ramp_buffer,npts);
S232('pushf',signal_ramp,npts);
S232('qpopf',TDT.ramp_buffer);

%% ZERO PLAY BUFFERS
zero_play_buffers;

%% START SEQUENCED PLAY
init_sequenced_play;


%% MAIN LOOP
seekpos=0;
while(seekpos < TDT.npts_total_play)
    
    % WAIT FOR LAST BUFFER TO FINISH
    while(check_play(TDT.nPlayChannels,[BUF_L1 BUF_R1])); end;
    
    tic;
    
    % SET FLAGS
    SignalPlayFlag=0;
    if(signalScale>0)
        readflag=1;
    elseif(readflag)
        readflag=0;
        SignalPlayFlag=1;
    end
    % COUNTDOWN (in seconds) TO NEXT TEST TRIAL
    cntdown=update_countdown(cnt,Signalcnt);
    % DISPLAY SESSION INFO
    disp_session_info(cntdown,seekpos);
    % PREPARE ADAPTOR
    if(PDR.flag_adapt)
        [adapt_state, adapt_left, adapt_right]=adaptor_filter(adapt_state);
    end
    % TEST TRIAL SCALE
    if(cnt==PDR.isi_buf)
        loc=PDR.TEST_loc_sequence(Signalcnt);
        if(loc~=0)
            signalScale=PDR.TEST_scale_sequence(Signalcnt);
        else
            signalScale=0;
        end
    end
    % LEFT CHANNEL BUFFER
    update_buffer(BUF_L1,adapt_left,test_left,cnt,Signalcnt,signalScale);
    % RIGHT CHANNEL BUFFER
    update_buffer(BUF_R1,adapt_right,test_right,cnt,Signalcnt,signalScale);
    % UPDATE ISI COUNTER AND SIGNAL COUNT
    if(cnt==PDR.isi_buf)
        cnt=round(TDT.max_signal_jitter.*rand);
        Signalcnt=Signalcnt+1;
    else
        cnt=cnt+1;
    end
    % RECORD PDR TRACE
    if(PDR.record)
        % First Record Channel:
        ch=1; buf=1;
        session.last_buffer=record_buffer(TDT.rec_buffers(ch,buf),TDT.dec_buffers(ch,buf),SignalPlayFlag,TDT.display_flag);
        session.test_flag=SignalPlayFlag;
        sessionPlots2('Update Trace Plot');
        % Second Record Channel:
        ch=2; buf=1;
        record_buffer(TDT.rec_buffers(ch,buf),TDT.dec_buffers(ch,buf),SignalPlayFlag,0);
    end
    % PROCESSING TIME
    if(TDT.disp_proc_time)
        t1=toc;
        figure(h1); delete(session.txt(4)); axis off;
        session.txt(4)= text(.01,.3,sprintf('Processing Time: %.3f seconds',t1),'FontSize',10);
        drawnow;
    end
        
    % UPDATE SEEK POSITION
    seekpos = seekpos + PDR.buf_pts;
    
    if(seekpos<TDT.npts_total_play)
        
        % WAIT FOR LAST BUFFER TO FINISH
        while(check_play(TDT.nPlayChannels,[BUF_L2 BUF_R2])); end;
        
        tic;
        
        % SET FLAGS
        SignalPlayFlag=0;
        if(signalScale>0)
            readflag=1;
        elseif(readflag)
            readflag=0;
            SignalPlayFlag=1;
        end
        % COUNTDOWN (in seconds) TO NEXT TEST TRIAL
        cntdown=update_countdown(cnt,Signalcnt);
        % DISPLAY SESSION INFO
        disp_session_info(cntdown,seekpos);
        % PREPARE ADAPTOR
        if(PDR.flag_adapt)
            [adapt_state, adapt_left, adapt_right]=adaptor_filter(adapt_state);
        end
        % TEST TRIAL SCALE
        if(cnt==PDR.isi_buf)
            loc=PDR.TEST_loc_sequence(Signalcnt);
            if(loc~=0)
                signalScale=PDR.TEST_scale_sequence(Signalcnt);
            else
                signalScale=0;
            end
        end
        % LEFT CHANNEL BUFFER
        update_buffer(BUF_L2,adapt_left,test_left,cnt,Signalcnt,signalScale);
        % RIGHT CHANNEL BUFFER
        update_buffer(BUF_R2,adapt_right,test_right,cnt,Signalcnt,signalScale);
        % RECORD PDR TRACE
        if(PDR.record)
            % First Record Channel:
            ch=1; buf=2;
            session.last_buffer=record_buffer(TDT.rec_buffers(ch,buf),TDT.dec_buffers(ch,buf),SignalPlayFlag,TDT.display_flag);
            session.test_flag=SignalPlayFlag;
            sessionPlots2('Update Trace Plot');
            % Second Record Channel:
            ch=2; buf=2;
            record_buffer(TDT.rec_buffers(ch,buf),TDT.dec_buffers(ch,buf),SignalPlayFlag,0);
        end
        % UPDATE ISI COUNTER AND SIGNAL COUNT
        if(cnt=PDR.isi_buf)
            cnt=round(TDT.max_signal_jitter.*rand);
            Signalcnt=Signalcnt+1;
        else
            cnt=cnt+1;
        end
        % UPDATE SEEK POSITION
        seekpos = seekpos + PDR.buf_pts;
        % CHECK IF CORRECT BUFFERS ARE PLAYING
        out=check_play(TDT.nPlayChannels,[BUF_L1 BUF_R1]);
        if(out==-1)
            disp(sprintf('Got %.2f percentof the way',seekpos/TDT.npts_total_play));
            disp('APcard too slow? or outFNs incorrect?');
            break;
        end
        % PROCESSING TIME
        if(TDT.disp_proc_time)
            t1=toc;
            figure(h1); delete(session.txt(4)); axis off;
            session.txt(4)= text(.01,.3,sprintf('Processing Time: %.3f seconds',t1),'FontSize',10);
            drawnow;
        end
    end
    
    if(Signalcnt>TDT.ntrials); Signalcnt=0; end;
    
end

%% WAIT FOR LAST BUFFERS TO FINISH
while(S232('playseg',TDT.din)==BUF_R2 || S232('playseg',TDT.din)==BUF_L2); end;
TDT_flush;

%% SUBROUTINES

function disp_session_info(cntdown,seekpos)
global TDT session
figure(session.hPlot);
subplot(session.hSessionInfo);
delete(session.txt(1)); delete(session.txt(2)); delete(session.txt(3)); axis off;
elapsed_time=seekpos*(TDT.srate/1e6);
min=floor(elapsed_time/60); sec=elapsed_time-(60*min);
session.txt(1)=session.txt(.01,.9,sprintf('ELAPSED TIME:    %i minutes   %.2f seconds',min,sec),'FontSize',12);
rem_time = TDT.npts_total_play*(SRATE/1E6) - elapsed_time;
min=floor(rem_time/60); sec=rem_time-(60*min);
session.txt(2)=session.txt(.01,.7,sprintf('REMAINING TIME:  %i minutes   %.2f seconds',min,sec),'FontSize',12);
min=floor(cntdown/60); sec=cntdown-(60*min);
session.txt(3)=session.txt(.01,.5,sprintf('NEXT TEST TRIAL: %i minutes   %.2f seconds',min,sec),'FontSize',12);
drawnow;

function cntdown=update_countdown(cnt,Signalcnt)
global PDR TDT
cntdown=(PDR.isi_buf-cnt)*(PDR.buf_pts*(TDT.srate/1e6));
if(PDR.TEST_scale_sequence(Signalcnt)==0)
    for(j=1:(PDR.ntrials-Signalcnt))
        cntdown=cntdown+(PDR.isi_buf+1)*(PDR.buf_pts*TDT.srate/1e6);
        if(PDR.TEST_scale_sequence(j+Signalcnt)>0)
            break;
        end
    end
end

function [adapt_state, adapt_left, adapt_right]=adaptor_filter(adapt_state)
global PDR HRTF
% get new set of pseudorandomly generated numbers and run through a FIR filter
rand('state',adapt_state);
new_buffer=rand(1,PDR.buf_pts);
adapt_state=rand('state');
[filtered_buffer, circ_buffer] = circ_fir(circ_buffer,new_buffer,PDR.ADAPT_coefs);
% filter adaptor with left/right HRTF coefficients
[adapt_left, left_buffer] = circ_fir(left_buffer,filtered_buffer,HRTF.AdaptL);
[adapt_right, right_buffer] = circ_fir(right_buffer,filtered_buffer,HRTF.AdaptR);

function update_buffer(BUF_ID,adaptor,test,cnt,Signalcnt,signalScale)
global PDR session
S232('dropall');
if(PDR.flag_adapt)
    S232('pushf',adaptor,PDR.buf_pts); scale(PDR.ADAPT_scale);
    if(cnt==PDR.isi_buf)
        S232('qpushf',PDR.ADAPT_ramp); S232('mult');
    end
else
    S232('dpush',PDR.buf_pts); S232('value',0);
end
if(cnt==PDR.isi_buf)
    loc=PDR.TEST_loc_sequence(Signalcnt);
    if(loc~=0)
        % plot a marker on trial sequence plot
        x=Signalcnt+1; y=loc;
        figure(session.hTrialPlot);
        if exist('hMark')
            delete(hMark);
        end
        hMark=plot(x,y,'MarkerSize',12,'Marker','s',...
            'MarkerFaceColor','none','MarkerEdgeColor','w');
        % push test onto stack and add to the ramped adaptor buffer:
        S232('qpushf',test);
        S232('scale',signalScale);
        S232('add'); % add to adaptor buffer
    end
end
S232('qpop16',BUF_ID);

function out=check_play(nPlayChannels,BUFFERS)
out=1;
for ch=1:nPlayChannels
    if(S232('playseg',ch)~=BUFFERS(ch))
        out=0; % error
    end
end

function out=check_params
global TDT PDR
out=1;
% check decimation factor
div=PDR.buf_pts/2^PDR.decimationfactor;
if(round(div)~=div)
    h = warndlg('buf_pts and decimationfactor incompatible');
    uiwait(h);
    out=-1;
end

% check npts total play
div=PDR.npts_totalplay/PDR.buf_pts;
if(round(div)~=div)
    h = warndlg('buf_pts and npts total play incompatible');
    uiwait(h);
    out=-1;
end

% check HRTF arrays
if(size(HRTF.AdaptL,1)~=1 && size(HRTF.AdaptL,2)~=PDR.HRTF_nlines)
    h = warndlg('Adaptor HRTF arrays must be 1 x HRTF_nlines');
    uiwait(h);
    out=-1;
end
if(size(HRTF.AdaptR,1)~=1 && size(HRTF.AdaptR,2)~=PDR.HRTF_nlines)
    h = warndlg('Adaptor HRTF arrays must be 1 x HRTF_nlines');
    uiwait(h);
    out=-1;
end
if(size(HRTF.TestL,1)~=1 && size(HRTF.TestL,2)~=PDR.HRTF_nlines)
    h = warndlg('Test HRTF arrays must be HRTF_nlines x # Test Locs');
    uiwait(h);
    out=-1;
end
if(size(HRTF.TestR,1)~=1 && size(HRTF.TestR,2)~=PDR.HRTF_nlines)
    h = warndlg('Test HRTF arrays must be HRTF_nlines x # Test Locs');
    uiwait(h);
    out=-1;
end

% check ramp vector
if(size(PDR.ADAPT_ramp,1)~=1 && size(PDR.ADAPT_ramp,2)~=PDR.buf_pts)
    h = warndlg('Ramp must be 1 x buf_pts');
    uiwait(h);
    out=-1;
end

% check test sound
if(size(PDR.TEST_sound,1)~=1 && size(PDR.TEST_sound,2)~=PDR.buf_pts)
    h = warndlg('Test sound must be 1 x buf_pts');
    uiwait(h);
    out=-1;
end

% check scale sequence
if(size(PDR.TEST_scale_sequence,1)~=1 && size(PDR.TEST_scale_sequence,2)<PDR.n_trials)
    h = warndlg('Scale sequence must be a 1 x n (n>=ntrials), row vector');
    uiwait(h);
    out=-1;
end

% check location sequence
if(size(PDR.TEST_loc_sequence,1)~=1 && size(PDR.TEST_loc_sequence,2)<PDR.n_trials)
    h = warndlg('Location sequence must be a 1 x n (n>=ntrials), row vector');
    uiwait(h);
    out=-1;
end


