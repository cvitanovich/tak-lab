function LDS_PDR_loop
% function with subroutines for double buffered LDS PDR expt
global TDT PDR session

% DEBUGGING FLAG TO SKIP TDT INSTRUCTIONS
DEBUG=PDR.DEBUG;
if(DEBUG)
    PDR.SOUNDS_speaker_scales_lag=30000*ones(1,length(PDR.location_seq));
    PDR.SOUNDS_speaker_scales_lead=PDR.SOUNDS_speaker_scales_lag;
    PDR.SOUNDS_lead_attens=zeros(1,length(PDR.location_seq));
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
TDT.outFN{1}=[PDR.filename '_REC1.vrt'];
TDT.outFN{2}=[PDR.filename '_REC2.vrt'];
TDT.ntrials=PDR.ntrials;
TDT.srate=1e6 / TDT.Fs;
TDT.display_flag=1; % flag to display trace during session
TDT.max_signal_jitter=PDR.TEST_trial_jitter;
TDT.disp_proc_time=1; % flag to display processing time for each buffer segment

% speaker ID of habituating lag location:
HAB_LOC=find(PDR.SOUNDS_azimuths==PDR.LAG_hab_pos);


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


if(~DEBUG) % for debugging without the TDT
    %% INITIALIZE TDT
    out=TDT_init;
    if(out==-1); return; end;
    %% INITIALIZE BUFFERS
    TDT=TDT_buffers(TDT);
    LEAD_PLAY = [TDT.stim_buffers{1}(1) TDT.stim_buffers{1}(2)];
    LAG_PLAY = [TDT.stim_buffers{2}(1) TDT.stim_buffers{2}(2)];
    REC_A = [TDT.rec_buffers{1}(1) TDT.rec_buffers{1}(2)];
    REC_B = [TDT.rec_buffers{2}(1) TDT.rec_buffers{2}(2)];
    DEC_A = [TDT.dec_buffers{1}(1) TDT.dec_buffers{1}(2)];
    DEC_B = [TDT.dec_buffers{2}(1) TDT.dec_buffers{2}(2)];
    %% INITIALIZE PD1
    PD1_init(TDT);
    %% INIT SIGNAL SWITCHERS
    S232('SS1clear',1); % LAG SS1
    S232('SS1mode',1,'QUAD_2_1'); % inputs 1, 3, 5, 7 => outputs A,B,C,D
    S232('SS1select',1,0,1); % Default lag output is A (hab location)
    % INIT ATTENS
    S232('PA4atten',1,0);
    S232('PA4atten',2,0);
    %% ZERO PLAY BUFFERS
    zero_play_buffers(TDT);
    %% START SEQUENCED PLAY
    init_sequenced_play(TDT);
end

%% MAIN LOOP
seekpos=0;
while(seekpos < TDT.npts_total_play)
    % WAIT FOR LAST BUFFER TO FINISH
    if(~DEBUG)
        while(check_play(TDT.nPlayChannels,[LEAD_PLAY(1) LAG_PLAY(1)])); end;
    end
    
    tic;
    % SET FLAGS
    SignalPlayFlag=0;
    if(signalScale>0);
        readflag=1;
    elseif(readflag>0)
        readflag=0;
        SignalPlayFlag=1;
    end
    
    % COUNTDOWN (in seconds) TO NEXT TEST TRIAL
    cntdown=update_countdown(cnt,Signalcnt);
    % DISPLAY SESSION INFO
    disp_session_info(cntdown,seekpos);
    
    % re-loading #1 playbuffers ... LEAD to chanA and LAG to chanB
    if(cnt==PDR.isi_buf) % stimulus trial
        % jitter trial presentation
        cnt=0;
        if(TDT.max_signal_jitter > 0)
            cnt=round(2*TDT.max_signal_jitter*rand-TDT.max_signal_jitter); % +/- max
        end
        % update location
        loc=PDR.location_seq(Signalcnt);
        
        % TDT instructions for this buffer:
        if(~DEBUG)
            S232('dropall');
            S232('PA4atten',1,PDR.SOUNDS_lead_attens(loc-1)); % lead channel
            S232('PA4atten',2,0); % lag ch
            
            % set SS1 units
            switch_speaker(loc);
            
            % select roved stimuli for buffer insertion
            S232('qpushf',PDR.LEAD_sounds{PDR.SOUNDS_rove_sequence(Signalcnt)});
            S232('scale',PDR.SOUNDS_speaker_scales_lead(loc-1));
            S232('qpop16',LEAD_PLAY(1));
            S232('qpushf',PDR.LAG_sounds{PDR.SOUNDS_rove_sequence(Signalcnt)});
            S232('scale',PDR.SOUNDS_speaker_scales_lag(loc-1));
            S232('qpop16',LAG_PLAY(1));
        end
        
        % plot a marker on trial sequence plot
        session.trialcnt=Signalcnt;
        session.trialval=PDR.SOUNDS_azimuths(loc);
        sessionPlots('Update Trial Plot');
        % increment signal count
        Signalcnt=Signalcnt+1;
        
    else	% no stimulus
        signalScale = 0;
        cnt=cnt+1; % increment ISI count
        if(~DEBUG)
            S232('dpush',PDR.buf_pts);
            S232('value',0);
            S232('qpop16',LEAD_PLAY(1));
            S232('dpush',PDR.buf_pts);
            S232('value',0);
            S232('qpop16',LEAD_PLAY(1));
        end
    end
    
    % RECORD PDR TRACE
    if(record && ~DEBUG)
        % First Record Channel:
        ch=1; buf=1;
        session.last_buffer=record_buffer(ch, REC_A(buf),DEC_A(buf),SignalPlayFlag,TDT.display_flag);
        if(SignalPlayFlag)
            if(PDR.location_seq(Signalcnt)==HAB_LOC)
                session.test_flag=1;
            else
                session.test_flag=Inf;
            end
        else
            session.test_flag=0;
        end
        sessionPlots('Update Trace Plot');
        % Second Record Channel:
        ch=2; buf=1;
        record_buffer(ch, REC_B(buf),DEC_B(buf),SignalPlayFlag,0);
        S232('dropall');
    end
    
    % UPDATE SEEK POSITION
    seekpos = seekpos + PDR.buf_pts;
    
    if(DEBUG)
        pause(PDR.buf_pts/PDR.stim_Fs);
    end
    % PROCESSING TIME
    if(TDT.disp_proc_time)
        t=toc;
        session.proc_time(end+1)=t;
        sessionPlots('Update Session Info');
    end
    
    if(seekpos<TDT.npts_total_play)
        
        % WAIT FOR LAST BUFFER TO FINISH
        if(~DEBUG)
            while(check_play(TDT.nPlayChannels,[LEAD_PLAY(2) LAG_PLAY(2)])); end;
        end
        
        tic;
        % SET FLAGS
        SignalPlayFlag=0;
        if(signalScale>0);
            readflag=1;
        elseif(readflag>0)
            readflag=0;
            SignalPlayFlag=1;
        end
        
        % COUNTDOWN (in seconds) TO NEXT TEST TRIAL
        cntdown=update_countdown(cnt,Signalcnt);
        % DISPLAY SESSION INFO
        disp_session_info(cntdown,seekpos);
        
        % re-loading #1 playbuffers ... LEAD to chanA and LAG to chanB
        if(cnt==PDR.isi_buf) % stimulus trial
            % jitter trial presentation
            cnt=0;
            if(TDT.max_signal_jitter > 0)
                cnt=round(2*TDT.max_signal_jitter*rand-TDT.max_signal_jitter); % +/- max
            end
            % update location
            loc=PDR.location_seq(Signalcnt);
            
            % TDT instructions for this buffer:
            if(~DEBUG)
                S232('dropall');
                S232('PA4atten',1,PDR.SOUNDS_lead_attens(loc-1)); % lead channel
                S232('PA4atten',2,0); % lag ch
                % set SS1 units
                switch_speaker(loc);
                
                % select roved stimuli for buffer insertion
                S232('qpushf',PDR.LEAD_sounds{PDR.SOUNDS_rove_sequence(Signalcnt)});
                S232('scale',PDR.SOUNDS_speaker_scales_lead(loc-1));
                S232('qpop16',BUF_A2);
                S232('qpushf',PDR.LAG_sounds{PDR.SOUNDS_rove_sequence(Signalcnt)});
                S232('scale',PDR.SOUNDS_speaker_scales_lag(loc-1));
                S232('qpop16',BUF_B2);
            end
            
            % plot a marker on trial sequence plot
            session.trialcnt=Signalcnt;
            session.trialval=PDR.SOUNDS_azimuths(loc);
            sessionPlots('Update Trial Plot');
            % increment signal count
            Signalcnt=Signalcnt+1;
            
        else	% no stimulus
            signalScale = 0; cnt=cnt+1;
            if(~DEBUG)
                S232('dpush',PDR.buf_pts);
                S232('value',0);
                S232('qpop16',BUF_A2);
                S232('dpush',PDR.buf_pts);
                S232('value',0);
                S232('qpop16',BUF_B2);
            end
        end
        
        % RECORD PDR TRACE
        if(record && ~DEBUG)
            % First Record Channel:
            ch=1; buf=2;
            session.last_buffer=record_buffer(ch, REC_A(buf),DEC_A(buf),SignalPlayFlag,TDT.display_flag);
            if(SignalPlayFlag)
                if(PDR.location_seq(Signalcnt)==HAB_LOC)
                    session.test_flag=1;
                else
                    session.test_flag=Inf;
                end
            else
                session.test_flag=0;
            end
            sessionPlots('Update Trace Plot');
            % Second Record Channel:
            ch=2; buf=2;
            record_buffer(ch, REC_B(buf),DEC_B(buf),SignalPlayFlag,0);
            S232('dropall');
        end
        
        % UPDATE SEEK POSITION
        seekpos = seekpos + PDR.buf_pts;
        
        % CHECK IF CORRECT BUFFERS ARE PLAYING
        if(~DEBUG)
            if(~check_play(TDT.nPlayChannels,[LEAD_PLAY(1) LAG_PLAY(1)]))
                disp(sprintf('Got %.2f percentof the way',seekpos/TDT.npts_total_play));
                disp('APcard too slow? or outFNs incorrect?');
                break;
            end
        end
        
        % PROCESSING TIME
        if(TDT.disp_proc_time)
            t=toc;
            session.proc_time(end+1)=t;
            sessionPlots('Update Session Info');
        end
        
        if(DEBUG)
            pause(PDR.buf_pts/PDR.stim_Fs);
        end
        
        if(Signalcnt>TDT.ntrials); break; end;
        
    end; % end second buffer segment
    
end

%% WAIT FOR LAST BUFFERS TO FINISH
if(~DEBUG)
    while(S232('playseg',TDT.din)==LAG_PLAY(2) || S232('playseg',TDT.din)==LEAD_PLAY(2)); end;
    TDT_flush;
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
sessionPlots('Update Session Info');

function cntdown=update_countdown(cnt,Signalcnt)
global PDR TDT
cntdown=(PDR.isi_buf-cnt)*(PDR.buf_pts*(TDT.srate/1e6));
for(j=1:(PDR.ntrials-Signalcnt))
    if(PDR.location_seq(j+Signalcnt)~=PDR.hab_loc_id)
        break;
    else
        cntdown=cntdown+(PDR.isi_buf+1)*(PDR.buf_pts*TDT.srate/1e6);
    end
end

function out=check_play(nPlayChannels,BUFFERS)
out=true;
for ch=1:nPlayChannels
    if(S232('playseg',ch)~=BUFFERS(ch))
        out=false; % error
    end
end

function switch_speaker(loc)
% switches speaker output of lag sound
S232('SS1clear',1); S232('SS1clear',2);
if(loc<6)
    ss_id=1; % use left SS1
    out_port=loc-2; % decrement by 2 for output selection
else
    ss_id=2; % right SS1
    out_port=loc-6; % decrement by 6 for output selection
end
S232('SS1mode',ss_id,'QUAD_2_1');
S232('SS1select',ss_id,out_port,1);

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