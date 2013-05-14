 function spatialPDR
% A FUNCTION TO RUN A SPATIAL PDR EXPERIMENT WITH AN ADAPTOR
global PDR TDT HRTF session


% DEBUGGING FLAG TO SKIP TDT INSTRUCTIONS
DEBUG=PDR.DEBUG;
if(DEBUG)
 	% SET ATTENS & SCALES, maybe??? (FIX)
end

%% INITIATE TDT PARAMETERS HERE
TDT.nPlayChannels=2;
TDT.playpts = {[PDR.buf_pts PDR.buf_pts],[PDR.buf_pts PDR.buf_pts]};
if(PDR.record>0)
    TDT.nRecChannels=2;
    TDT.recpts=TDT.playpts;
else
    TDT.nRecChannels=0;
end
TDT.dec_factor=PDR.decimationfactor;
TDT.din = 1;
TDT.Fs = PDR.stim_Fs;
TDT.npts_total_play=PDR.npts_totalplay;
% data storage files:
if PDR.virtual
    type='.vrt';
else
    type='.frf';
end
TDT.outFN{1}=[PDR.filename '_REC1' type];
TDT.outFN{2}=[PDR.filename '_REC2' type];
TDT.ntrials=PDR.ntrials;
TDT.srate=1e6 / TDT.Fs;
TDT.display_flag=1; % flag to display trace during session
TDT.max_signal_jitter=PDR.jitter;
TDT.disp_proc_time=1; % flag to display processing time for each buffer segment
TDT.hab_id=PDR.hab_id;
% test stimulus variable:
TDT.stim_left=zeros(1,PDR.buf_pts);
TDT.stim_right=TDT.stim_left;

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


%% ADAPTOR FILTERING PARAMS
% recording adaptor sequence:
PDR.ADAPT_state_seq=[];
% initial state (picked randomly from list of adaptor states)
rand_idx=round((length(PDR.ADAPT_state_list)-1)*rand+1);
adapt_state = PDR.ADAPT_state_list(rand_idx);
PDR.ADAPT_state_seq(end+1)=rand_idx;
% circular buffer for continuous filtered adaptor:
CIRC_BUFS.adaptor=zeros(1,(length(PDR.ADAPT_coefs)+PDR.buf_pts));
% circular buffers for HRTF filtering (left/right):
CIRC_BUFS.left=CIRC_BUFS.adaptor; CIRC_BUFS.right=CIRC_BUFS.adaptor;

%% FILTER TEST SOUNDS WITH HRTFS
% buffer assignments for TDT
% zeros added later to get buffer size = PDR.buf_pts!!!!
TESTLEFT=zeros(PDR.TEST_nlocs,PDR.TEST_stim_pts+2*PDR.HRTF_nTaps);
TESTRIGHT=TESTLEFT;
filtered_test_left=zeros(1,PDR.TEST_stim_pts+2*PDR.HRTF_nTaps);
filtered_test_right=filtered_test_left;
% loop through test sounds
for(j=1:PDR.TEST_nlocs)
    % filter each test sound with HRTFS & Store on AP2 Card
    filtered_test_left=filter(HRTF.TestL(:,j),1,PDR.TEST_sound);
    filtered_test_right=filter(HRTF.TestR(:,j),1,PDR.TEST_sound);
    % get rms of stimulus (ignoring empty padding for HRIR filtering)
    tmp_Lt=find(filtered_test_left~=0); tmp_Rt=find(filtered_test_right~=0);
    rms_val=(sqrt(mean(filtered_test_left(tmp_Lt(1):tmp_Lt(end)).^2))+...
        sqrt(mean(filtered_test_right(tmp_Rt(1):tmp_Rt(end)).^2)))/2;
    filtered_test_left = (PDR.TEST_target_rms/rms_val).*filtered_test_left;
    filtered_test_right = (PDR.TEST_target_rms/rms_val).*filtered_test_right;
    TESTLEFT(j,:)=filtered_test_left;
    TESTRIGHT(j,:)=filtered_test_right;
end
clear filtered_test_left filtered_test_right

if(~DEBUG) % for debugging without the TDT
    %% INITIALIZE TDT
    out=TDT_init;
    if(out==-1); return; end;
    
    %% INITIALIZE BUFFERS
    TDT=TDT_buffers(TDT);
    LEFT_PLAY = [TDT.stim_buffers{1}(1) TDT.stim_buffers{1}(2)];
    RIGHT_PLAY = [TDT.stim_buffers{2}(1) TDT.stim_buffers{2}(2)];
    REC_A = [TDT.rec_buffers{1}(1) TDT.rec_buffers{1}(2)];
    REC_B = [TDT.rec_buffers{2}(1) TDT.rec_buffers{2}(2)];
    DEC_A = [TDT.dec_buffers{1}(1) TDT.dec_buffers{1}(2)];
    DEC_B = [TDT.dec_buffers{2}(1) TDT.dec_buffers{2}(2)];
    
    %% INITIALIZE PD1
    PD1_init(TDT);
    
    
    %% SET ATTENS
    TDT.attens=[PDR.base_atten PDR.base_atten];
    TDT_attens(TDT);
    
    %% SIGNAL RAMP BUFFER
    TDT.ramp_buffer=TDT.n_total_buffers+1;
    S232('allotf',TDT.ramp_buffer,PDR.buf_pts);
    S232('pushf',PDR.ADAPT_ramp,PDR.buf_pts);
    S232('qpopf',TDT.ramp_buffer);
    TDT.n_total_buffers=TDT.n_total_buffers+1;
    
    %% ZERO PLAY BUFFERS
    zero_play_buffers(TDT);
    
    %% START SEQUENCED PLAY
    init_sequenced_play(TDT);
    
    
end

%% MAIN LOOP
seekpos=0;
while(seekpos < TDT.npts_total_play)
    if(~DEBUG)
        % WAIT FOR LAST BUFFER TO FINISH
        while(check_play(TDT.nPlayChannels,[LEFT_PLAY(1) RIGHT_PLAY(1)])); end;
    end
    
    tic;
    
    % SET FLAGS
    SignalPlayFlag=0;
    if(signalScale>0)
        readflag=1;
    elseif(readflag>0)
        readflag=0;
        SignalPlayFlag=1;
    end
    % COUNTDOWN (in seconds) TO NEXT TEST TRIAL
    cntdown=update_countdown(cnt,Signalcnt);
    % DISPLAY SESSION INFO
    disp_session_info(cntdown,seekpos);
    
    % TEST TRIAL SCALE
    
    TDT.stim_left=zeros(1,PDR.buf_pts);
    TDT.stim_right=zeros(1,PDR.buf_pts);
    session.stim_left=zeros(1,session.stim_pts);
    session.stim_right=zeros(1,session.stim_pts);
    if(cnt==PDR.isi_buf)
        loc=PDR.TEST_loc_sequence(Signalcnt);
        if(loc~=0)
            signalScale=PDR.TEST_scale_sequence(Signalcnt);
            TDT.stim_left(PDR.TEST_start_pt:PDR.TEST_stop_pt)=TESTLEFT(loc,:);
            TDT.stim_right(PDR.TEST_start_pt:PDR.TEST_stop_pt)=TESTRIGHT(loc,:);
        else % not playing test sound in this trial!
            signalScale=0;
        end
        % scale test
        TDT.stim_left=round(signalScale.*TDT.stim_left);
        TDT.stim_right=round(signalScale.*TDT.stim_right);
        % plot a marker on trial sequence plot
        session.trialcnt=Signalcnt;
        session.trialval=loc;
        sessionPlots3('Update Trial Plot');
        % plot sound buffers:
        session.stim_left=TDT.stim_left(PDR.TEST_start_pt:PDR.TEST_stop_pt);
        session.stim_right=TDT.stim_right(PDR.TEST_start_pt:PDR.TEST_stop_pt);
        sessionPlots3('Update Stim Plot');
    end
    
    
    % SETUP ADAPTOR
    
    if(cnt~=PDR.isi_buf && PDR.flag_adapt>0)
        [adapt_state, adapt_left, adapt_right, CIRC_BUFS]=adaptor_filter(adapt_state,CIRC_BUFS);
    elseif(cnt==PDR.isi_buf && PDR.flag_adapt>0)
        % clear circular buffers & set new seed value
        rand_idx=round((length(PDR.ADAPT_state_list)-1)*rand+1);
        adapt_state = PDR.ADAPT_state_list(rand_idx);
        PDR.ADAPT_state_seq(end+1)=rand_idx; % save state sequence
        CIRC_BUFS.adaptor=zeros(1,(length(PDR.ADAPT_coefs)+PDR.buf_pts));
        CIRC_BUFS.left=CIRC_BUFS.adaptor; CIRC_BUFS.right=CIRC_BUFS.adaptor;
        % paste prior/new buffers together (the trial ramp will remove any
        % possible discontinuities in the adaptor)
        [adapt_state, aleft_new, aright_new, CIRC_BUFS]=adaptor_filter(adapt_state,CIRC_BUFS);
        adapt_left=[adapt_left(1:floor(end/2)) ...
            aleft_new(ceil(end/2):end)];
        adapt_right=[adapt_left(1:floor(end/2)) ...
            aright_new(ceil(end/2):end)];
    end
    
    if(~DEBUG)
        % LEFT CHANNEL BUFFER
        update_buffer(LEFT_PLAY(1),adapt_left,TDT.stim_left,cnt,Signalcnt,signalScale);
        % RIGHT CHANNEL BUFFER
        update_buffer(RIGHT_PLAY(1),adapt_right,TDT.stim_right,cnt,Signalcnt,signalScale);
    end
    
    % UPDATE ISI COUNTER AND SIGNAL COUNT
    if(cnt==PDR.isi_buf)
        cnt=round(2*TDT.max_signal_jitter*rand-TDT.max_signal_jitter); % +/- max
        Signalcnt=Signalcnt+1;
    else
        cnt=cnt+1;
    end
    % RECORD PDR TRACE
    if(record && ~DEBUG)
        % First Record Channel:
        ch=1; buf=1;
        session.last_buffer=record_buffer(ch, REC_A(buf),DEC_A(buf),SignalPlayFlag,TDT.display_flag);
        if(SignalPlayFlag==1)
            if(PDR.TEST_loc_sequence(Signalcnt)==TDT.hab_id)
                session.test_flag=1;
            else
                session.test_flag=Inf;
            end
        else
            session.test_flag=0;
        end
        sessionPlots3('Update Trace Plot');
        % Second Record Channel:
        ch=2; buf=1;
        record_buffer(ch, REC_B(buf),DEC_B(buf),SignalPlayFlag,0);
    end
    
    % UPDATE SEEK POSITION
    seekpos = seekpos + PDR.buf_pts;
    

    
    % CHECK FOR HALT
    if(session.HALT==1 && session.confirm_halt==0)
        confirmdlg;
    end
    if(session.HALT==1 && session.confirm_halt==1)
        seekpos=Inf;
    end
    

    % PROCESSING TIME
    if(TDT.disp_proc_time)
        t=toc;
        session.proc_time(end+1)=t;
        sessionPlots3('Update Session Info');
    end
    
    if(seekpos<TDT.npts_total_play)
        if(~DEBUG)
            % WAIT FOR LAST BUFFER TO FINISH
            while(check_play(TDT.nPlayChannels,[LEFT_PLAY(2) RIGHT_PLAY(2)])); end;
        end
        
        tic;
        
        % SET FLAGS
        SignalPlayFlag=0;
        if(signalScale>0)
            readflag=1;
        elseif(readflag>0)
            readflag=0;
            SignalPlayFlag=1;
        end
        % COUNTDOWN (in seconds) TO NEXT TEST TRIAL
        cntdown=update_countdown(cnt,Signalcnt);
        % DISPLAY SESSION INFO
        disp_session_info(cntdown,seekpos);
    
        
        % TEST TRIAL SCALE
        TDT.stim_left=zeros(1,PDR.buf_pts);
        TDT.stim_right=zeros(1,PDR.buf_pts);
        session.stim_left=zeros(1,session.stim_pts);
        session.stim_right=zeros(1,session.stim_pts);
        if(cnt==PDR.isi_buf)
            loc=PDR.TEST_loc_sequence(Signalcnt);
            if(loc~=0)
                signalScale=PDR.TEST_scale_sequence(Signalcnt);
                TDT.stim_left(PDR.TEST_start_pt:PDR.TEST_stop_pt)=TESTLEFT(loc,:);
                TDT.stim_right(PDR.TEST_start_pt:PDR.TEST_stop_pt)=TESTRIGHT(loc,:);
            else % not playing test sound in this trial!
                signalScale=0;
            end
            % scale test
            TDT.stim_left=round(signalScale.*TDT.stim_left);
            TDT.stim_right=round(signalScale.*TDT.stim_right);
            % plot a marker on trial sequence plot
            session.trialcnt=Signalcnt;
            session.trialval=loc;
            sessionPlots3('Update Trial Plot');
            % plot sound buffers:
            session.stim_left=TDT.stim_left(PDR.TEST_start_pt:PDR.TEST_stop_pt);
            session.stim_right=TDT.stim_right(PDR.TEST_start_pt:PDR.TEST_stop_pt);
            sessionPlots3('Update Stim Plot');
        end
        % SETUP ADAPTOR
        
        if(cnt~=PDR.isi_buf && PDR.flag_adapt>0)
            [adapt_state, adapt_left, adapt_right, CIRC_BUFS]=adaptor_filter(adapt_state,CIRC_BUFS);
        elseif(cnt==PDR.isi_buf && PDR.flag_adapt>0)
            % initialize circular buffers & set new seed value
            rand_idx=round((length(PDR.ADAPT_state_list)-1)*rand+1);
            adapt_state = PDR.ADAPT_state_list(rand_idx);
            PDR.ADAPT_state_seq(end+1)=rand_idx; % save state sequence
            CIRC_BUFS.adaptor=zeros(1,(length(PDR.ADAPT_coefs)+PDR.buf_pts));
            CIRC_BUFS.left=CIRC_BUFS.adaptor; CIRC_BUFS.right=CIRC_BUFS.adaptor;
            % paste prior/new buffers together (the trial ramp will remove any
            % possible discontinuities in the adaptor)
            [adapt_state, aleft_new, aright_new, CIRC_BUFS]=adaptor_filter(adapt_state,CIRC_BUFS);
            adapt_left=[adapt_left(1:floor(end/2)) ...
                aleft_new(ceil(end/2):end)];
            adapt_right=[adapt_left(1:floor(end/2)) ...
                aright_new(ceil(end/2):end)];
        end
        
        
        if(~DEBUG)
            % LEFT CHANNEL BUFFER
            update_buffer(LEFT_PLAY(2),adapt_left,TDT.stim_left,cnt,Signalcnt,signalScale);
            % RIGHT CHANNEL BUFFER
            update_buffer(RIGHT_PLAY(2),adapt_right,TDT.stim_right,cnt,Signalcnt,signalScale);
        end
        % UPDATE ISI COUNTER AND SIGNAL COUNT
        if(cnt==PDR.isi_buf)
            cnt=round(2*TDT.max_signal_jitter*rand-TDT.max_signal_jitter); % +/- max
            Signalcnt=Signalcnt+1;
        else
            cnt=cnt+1;
        end
        % RECORD PDR TRACE
        if(record && ~DEBUG)
            % First Record Channel:
            ch=1; buf=2;
            session.last_buffer=record_buffer(ch, REC_A(buf),DEC_A(buf),SignalPlayFlag,TDT.display_flag);
            if(SignalPlayFlag==1)
                if(PDR.TEST_loc_sequence(Signalcnt)==TDT.hab_id)
                    session.test_flag=1;
                else
                    session.test_flag=Inf;
                end
            else
                session.test_flag=0;
            end
            sessionPlots3('Update Trace Plot');
            % Second Record Channel:
            ch=2; buf=2;
            record_buffer(ch, REC_B(buf),DEC_B(buf),SignalPlayFlag,0);
        end
        
        
        % CHECK IF CORRECT BUFFERS ARE PLAYING
        if(~DEBUG)
            if(~check_play(TDT.nPlayChannels,[LEFT_PLAY(1) RIGHT_PLAY(1)]))
                str=[];
                for(k=1:TDT.nPlayChannels)
                    p=S232('playseg',k);
                    str=[str 'Ch ' num2str(k) ' playing buffer ' num2str(p) ' '];
                end
                disp(str);
                disp(sprintf('Got %.2f percent of the way',seekpos/TDT.npts_total_play));
                disp('APcard too slow? Or outFNs incorrect?');
                break;
            end
        end
        
        % CHECK FOR HALT
        if(session.HALT==1 && session.confirm_halt==0)
            confirmdlg;
        end
        if(session.HALT==1 && session.confirm_halt==1)
            seekpos=Inf;
        end
        
        % PROCESSING TIME
        if(TDT.disp_proc_time)
            t=toc;
            session.proc_time(end+1)=t;
            sessionPlots3('Update Session Info');
        end
    end
    
    if(Signalcnt>TDT.ntrials); break; end;
    
end

%% WAIT FOR LAST BUFFERS TO FINISH
if(~DEBUG)
    while(S232('playseg',TDT.din)==LEFT_PLAY(2) || S232('playseg',TDT.din)==RIGHT_PLAY(2)); end;
    TDT_flush(TDT);
end

if(session.HALT==1)
    disp('HALTING EARLY!!!');
    session.hHalt=warndlg('Halting Early!!!');
    uiwait(session.hHalt)
end


%% SUBROUTINES

function disp_session_info(cntdown,seekpos)
global TDT session
elapsed_time=seekpos*(TDT.srate/1e6);
min=floor(elapsed_time/60); sec=elapsed_time-(60*min);
session.elapsed_time=[min sec];
rem_time=TDT.npts_total_play*(TDT.srate/1E6) - elapsed_time;
min=floor(rem_time/60); sec=rem_time-(60*min);
session.rem_time=[min sec];
min=floor(cntdown/60); sec=cntdown-(60*min);
session.next_test_trial=[min sec];
sessionPlots3('Update Session Info');

function confirmdlg
global session
% dialog box to confirm a choice
session.Confirm_txt=uicontrol('Tag','confirm_txt','Style', 'text','String','Sure???',...
    'Units','normalized','FontSize',14,'Position',[0.02 0.5 0.05 0.05],...
    'BackgroundColor','w','ForegroundColor','k');
% buttons
figure(session.hFig); hold on;
yes='session.confirm_halt=1; session.HALT=1;';
nope='session.confirm_halt=0; session.HALT=0;';
nope=[nope ' h=findobj(''Tag'',''confirm''); h2=findobj(''Tag'',''confirm_txt'');'];
nope=[nope ' delete(h); delete(h2);'];
nope=[nope ' set(session.HALT_btn,''String'',''HALT'');'];
nope=[nope ' set(session.HALT_btn,''BackgroundColor'',''r'');'];
nope=[nope ' set(session.HALT_btn,''ForegroundColor'',''y'');'];
nope=[nope ' set(session.HALT_btn,''Callback'',''session.HALT=1;''); drawnow;'];
set(session.HALT_btn,'String','NO!',...
    'Units','normalized','FontSize',14,'Position',[0.02 0.4 0.05 0.05],...
    'BackgroundColor','w','ForegroundColor','k',...
    'Callback', nope);
session.Confirm_btn=uicontrol('Tag','confirm','Style', 'pushbutton','String','YES!',...
    'Units','normalized','FontSize',14,'Position',[0.02 0.3 0.05 0.05],...
    'BackgroundColor','w','ForegroundColor','k',...
    'Callback', yes);

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

function [adapt_state, adapt_left, adapt_right, CIRC_BUFS]=adaptor_filter(adapt_state,CIRC_BUFS)
global PDR HRTF
% get new set of pseudorandomly generated numbers and run through a FIR filter
rand('state',adapt_state);
new_buffer=rand(1,PDR.buf_pts);
adapt_state=rand('state');

[filtered_buffer, CIRC_BUFS.adaptor] = circ_fir(CIRC_BUFS.adaptor,new_buffer,PDR.ADAPT_coefs);
% filter adaptor with left/right HRTF coefficients
[adapt_left, CIRC_BUFS.left] = circ_fir(CIRC_BUFS.left,filtered_buffer,HRTF.AdaptL);
[adapt_right, CIRC_BUFS.right] = circ_fir(CIRC_BUFS.right,filtered_buffer,HRTF.AdaptR);

function update_buffer(BUF_ID,adaptor,test,cnt,Signalcnt,signalScale)
global PDR
S232('dropall');
if(PDR.flag_adapt)
    S232('pushf',adaptor,PDR.buf_pts); S232('scale',PDR.ADAPT_scale);
    if(cnt==PDR.isi_buf)
        S232('pushf',PDR.ADAPT_ramp,PDR.buf_pts); S232('mult');
    end
else
    S232('dpush',PDR.buf_pts); S232('value',0);
end
if(cnt==PDR.isi_buf)
    loc=PDR.TEST_loc_sequence(Signalcnt);
    if(loc~=0)
        % push test onto stack and add to the ramped adaptor buffer:
        S232('push16',test,PDR.buf_pts); % already scaled
        S232('add'); % add to adaptor buffer
    end
end
S232('qpop16',BUF_ID);

function out=check_play(nPlayChannels,BUFFERS)
out=true;
for ch=1:nPlayChannels
    if(S232('playseg',ch)~=BUFFERS(ch))
        out=false; % error
    end
end

function out=check_params
global TDT PDR HRTF
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

if(size(HRTF.TestL,2)~=PDR.TEST_nlocs && size(HRTF.TestL,1)~=PDR.HRTF_nlines)
    h = warndlg('Test HRTF arrays must be HRTF_nlines x # Test Locs');
    uiwait(h);
    out=-1;
end
if(size(HRTF.TestR,2)~=PDR.TEST_nlocs && size(HRTF.TestR,1)~=PDR.HRTF_nlines)
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


