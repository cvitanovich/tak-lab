function spatialPDR
% A FUNCTION TO RUN A SPATIAL PDR EXPERIMENT WITH AN ADAPTOR
global PDR TDT HRTF session

% DEBUGGING FLAG TO SKIP TDT INSTRUCTIONS
DEBUG=PDR.DEBUG;

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
if(PDR.virtual==1)
    type='.vrt';
else
    type='.frf';
end
TDT.outFN{1}=[PDR.data_path PDR.filename '_REC1' type];
TDT.outFN{2}=[PDR.data_path PDR.filename '_REC2' type];
TDT.ntrials=PDR.ntrials;
TDT.srate=1e6 / TDT.Fs; % srate is in usecs
TDT.display_flag=1; % flag to display trace during session
TDT.max_signal_jitter=PDR.jitter;
TDT.disp_proc_time=1; % flag to display processing time for each buffer segment
TDT.hab_id=PDR.hab_id; % ID for habituating trials (in this case... scale=zero)
% test stimulus variable:
TDT.stim_left=zeros(1,PDR.buf_pts);
TDT.stim_right=TDT.stim_left;
% buffer counter
TDT.buf_cnt=0;
TDT.total_buffers=TDT.npts_total_play*PDR.buf_pts;

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
    % important: since I set this up to use PD1 for HRTF filtering,
    % the play channels now correspond with (a) TEST and (b) ADAPTOR
    TEST_PLAY = [TDT.stim_buffers{1}(1) TDT.stim_buffers{1}(2)];
    ADAPT_PLAY = [TDT.stim_buffers{2}(1) TDT.stim_buffers{2}(2)];
    REC_A = [TDT.rec_buffers{1}(1) TDT.rec_buffers{1}(2)];
    REC_B = [TDT.rec_buffers{2}(1) TDT.rec_buffers{2}(2)];
    DEC_A = [TDT.dec_buffers{1}(1) TDT.dec_buffers{1}(2)];
    DEC_B = [TDT.dec_buffers{2}(1) TDT.dec_buffers{2}(2)];
    
    %% INITIALIZE PD1
    PD1_Route1(TDT); % routing commands for PD1

    %% LOAD HRTF COEFS TO PD1
    load_HRTF_coefs(HRTF);
    
    %% SET ATTENS
    TDT.attens=[PDR.base_atten PDR.base_atten];
    TDT_attens(TDT);
    
    if(PDR.flag_adapt>0)
        %% SIGNAL RAMP BUFFER
        TDT.ramp_buffer=TDT.n_total_buffers+1;
        S232('allotf',TDT.ramp_buffer,PDR.buf_pts);
        S232('pushf',PDR.ADAPT_ramp,PDR.buf_pts);
        S232('qpopf',TDT.ramp_buffer);
        TDT.n_total_buffers=TDT.n_total_buffers+1;
    end
end
   
%% store test sound on AP2
s232('dropall');
TDT.test_buf=TDT.n_total_buffers+1;
s232('allotf',TDT.test_buf,PDR.buf_pts);
s232('pushf',PDR.TEST_sound,PDR.buf_pts);
s232('qpopf',TDT.test_buf);
TDT.n_total_buffers=TDT.n_total_buffers+1;

%% store adaptor on AP2
s232('dropall');
TDT.adapt_buf=TDT.n_total_buffers+1;
adapt_pts=length(PDR.ADAPT_sound);
s232('allotf',TDT.adapt_buf,adapt_pts);
s232('pushf',PDR.ADAPT_sound,adapt_pts);
s232('qpopf',TDT.adapt_buf);
TDT.n_total_buffers=TDT.n_total_buffers+1;
% clear adaptor pos. spec
clear SPOS;
s232('dropall');

%% LINEAR GATE FOR FIRST/LAST ADAPTOR BUFFERS
% not to be confused with the ramp for the empty segment of each trial
TDT.gate_dur=10.0; % 10 ms gate duration
TDT.gate_pts=1E3*TDT.gate_dur/(TDT.srate); % gate duration in ms
step=1.0/TDT.gate_pts;
% create and store gates ahead of time
TDT.gate_on=TDT.n_total_buffers+1;
TDT.gate_off=TDT.n_total_buffers+2;
TDT.n_total_buffers=TDT.n_total_buffers+2;
s232('dpush',TDT.gate_pts);
s232('fill',0.0,step); % linear ramp for ON gate (0...1)
s232('allotf',TDT.gate_on,TDT.gate_pts);
s232('qpopf',TDT.gate_on);
s232('dpush',TDT.gate_pts);
s232('fill',(1.0-step),-step); % linear ramp for OFF gate (1...0)
s232('allotf',TDT.gate_off,TDT.gate_pts);
s232('qpopf',TDT.gate_off);
s232('dropall');

if(~DEBUG)
    %% ZERO PLAY BUFFERS
    zero_play_buffers(TDT);
    %% START SEQUENCED PLAY
    init_sequenced_play(TDT);
end

%% MAIN LOOP
seekpos=0; writeflag=0; readflag=0;
cycle=2; % cycle through double buffering seq
while(seekpos < TDT.npts_total_play)
 
    % cycle through double buffering sequence
    last=cycle;
    cycle=mod(cycle,2)+1;
    
    if(~DEBUG)
        % WAIT FOR LAST BUFFER TO FINISH
        while(check_play(TDT.nPlayChannels,[TEST_PLAY(cycle) ADAPT_PLAY(cycle)])); end;
    end
    
    % proc time
    tic;
    
    % update buffer counter
    TDT.buf_cnt=TDT.buf_cnt+1;
    
    
    % SET FLAGS
    SignalPlayFlag=0;
    if(writeflag==1) % currently playing trial buffer
        readflag=1; % set read flag
    elseif(readflag==1) % last buffer played was a trial buffer
        SignalPlayFlag=1; % insert a one at beginning of record buffer to indicate trial occured
        readflag=0; % reset read flag
    end
    
    % COUNTDOWN (in seconds) TO NEXT TEST TRIAL
    cntdown=update_countdown(cnt,Signalcnt);
    % DISPLAY SESSION INFO
    disp_session_info(cntdown,seekpos);
    
    % TEST TRIAL SCALE
    if(cnt==PDR.isi_buf)
        writeflag=1; % next buffer played will be a trial buffer
        signalScale=PDR.TEST_scale_sequence(Signalcnt);
        % plot a marker on trial sequence plot
        session.trialcnt=Signalcnt;
        sessionPlots_v4('Update Trial Plot');
    else
        signalScale=0;
        writeflag=0; % not writing trial buffer
    end

    % SETUP ADAPTOR
    if( ~exist('SPOS') || ( (SPOS+PDR.buf_pts+1) > adapt_pts ) )
        SPOS=0; % initialize specified portion of adaptor to use current buffer to position zero
    else
        SPOS=SPOS+PDR.buf_pts; % increment SPOS for next adaptor segment
    end

    % UPDATE PLAY BUFFERS
    if(~DEBUG)
        update_buffers(TEST_PLAY(cycle),ADAPT_PLAY(cycle),cnt,Signalcnt,signalScale,SPOS);
    end
  
    % RECORD PDR TRACE
    if(record && ~DEBUG)
        % First Record Channel:
        session.last_buffer=record_buffer(1, REC_A(cycle),DEC_A(cycle),TDT,SignalPlayFlag,TDT.display_flag);
        if(SignalPlayFlag==1)
            session.trial_flag=Inf; % no test sound (DEFAULT value)
            if((Signalcnt-1)>0)
                trial=Signalcnt-1;
                if(PDR.TEST_scale_sequence(trial)>0)
                    session.trial_flag=1; % test sound affirmative
                end
            end
        else
            session.trial_flag=0;
        end
        sessionPlots_v4('Update Trace Plot');
        % Second Record Channel:
        record_buffer(2, REC_B(cycle),DEC_B(cycle),TDT,SignalPlayFlag,0);
    end
    
    % UPDATE ISI COUNTER AND SIGNAL COUNT
    if(cnt==PDR.isi_buf)
        cnt=round(2*TDT.max_signal_jitter*rand-TDT.max_signal_jitter); % +/- max
        Signalcnt=Signalcnt+1;
    else
        cnt=cnt+1;
    end
    
    % CHECK IF CORRECT BUFFERS ARE PLAYING
    if(~DEBUG)
        if(~check_play(TDT.nPlayChannels,[TEST_PLAY(last) ADAPT_PLAY(last)]))
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
        sessionPlots_v4('Update Session Info');
    end
    if(Signalcnt>TDT.ntrials); break; end;
    
end

%% WAIT FOR LAST BUFFERS TO FINISH
if(~DEBUG)
    while(S232('playseg',TDT.din)==TEST_PLAY(last) || S232('playseg',TDT.din)==ADAPT_PLAY(last)); end;
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
sessionPlots_v4('Update Session Info');

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

function update_buffers(test_play,adapt_play,cnt,Signalcnt,signalScale,SPOS)
global PDR TDT

% prep adaptor buffer
S232('dropall');
if(PDR.flag_adapt>0)
    s232('qpushpartf',TDT.adapt_buf,SPOS,PDR.buf_pts); % pop correct adaptor segment to stack
    if(cnt==PDR.isi_buf)
        S232('qpushf',TDT.ramp_buffer); S232('mult'); % multiply adaptor with ramp if isi reached
    end
    if(TDT.buf_cnt==1)
        s232('qpushf',TDT.gate_on); % ramp ON first buffer in expt
        s232('mult');
    end
    if(TDT.buf_cnt==TDT.total_buffers)
        n=s232('topsize');
        s232('block',(n-TDT.gate_pts),(n-1));
        s232('qpushf',TDT.gate_off); % ramp OFF last buffer in expt
        s232('mult');
        s232('noblock');
    end
    S232('scale',PDR.ADAPT_scale); % scale adaptor
else
    S232('dpush',PDR.buf_pts); S232('value',0);
end
S232('qpop16',adapt_play);

% prep play buffer
S232('dropall');
if(cnt==PDR.isi_buf && signalScale~=0)
    % push test onto stack, SCALE and add to the ramped adaptor buffer:
    S232('qpushf',TDT.test_buf);
    S232('scale',round(signalScale));
else
    s232('dpush',PDR.buf_pts); S232('value',0);
end
S232('qpop16',test_play);

function out=check_play(nPlayChannels,BUFFERS)
out=true;
for ch=1:nPlayChannels
    if(S232('playseg',ch)~=BUFFERS(ch))
        out=false; % error
    end
end

function PD1_Route1(TDT)
% commands to set up routing in the TDT
SRATE = 1e6 / TDT.Fs;
S232('PD1clear',TDT.din);
S232('PD1fixbug',TDT.din);
S232('PD1nstrms',TDT.din,2,2);
S232('PD1srate',TDT.din,SRATE);
S232('PD1npts',TDT.din,TDT.npts_total_play);

S232('PD1clrsched',TDT.din);

% left earphone routing
src=[s232('DSPout',0) s232('DSPout',2)]; sf=[1 1];
s232('PD1addmult',TDT.din,src,sf,2,s232('DAC',0));

% right earphone routing
src=[s232('DSPout',1) s232('DSPout',3)]; sf=[1 1];
s232('PD1addmult',TDT.din,src,sf,2,s232('DAC',1));

% IREGS
S232('PD1addsimp',TDT.din,s232('IREG',0),s232('DSPin',0));
S232('PD1addsimp',TDT.din,s232('IREG',0),s232('DSPin',1));
S232('PD1addsimp',TDT.din,s232('IREG',1),s232('DSPin',2));
S232('PD1addsimp',TDT.din,s232('IREG',1),s232('DSPin',3));

% test sound routing
s232('PD1specIB',TDT.din,s232('IB',0),s232('IREG',0));

% adaptor routing
s232('PD1specIB',TDT.din,s232('IB',1),s232('IREG',1));

% pupillometer recording
s232('PD1specOB',TDT.din,s232('OB',0),s232('ADC',0));
s232('PD1specOB',TDT.din,s232('OB',1),s232('ADC',1));

% LED thresholds
S232('PD1setIO',TDT.din,0.01,9.99,0.01,9.99);

function load_HRTF_coefs(HRTF)
global PDR TDT
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

% load test coefs to DSP
s232('pushf',HRTF.TestL,PDR.HRTF_nlines);
s232('PreLoadRaw',TDT.din,s232('DSPid',0),MONO,STACK,[],[],1,1,1);
s232('pushf',HRTF.TestR,PDR.HRTF_nlines);
s232('PreLoadRaw',TDT.din,s232('DSPid',1),MONO,STACK,[],[],1,1,1);
% load adapt coefs to DSP
s232('pushf',HRTF.AdaptL,PDR.HRTF_nlines);
s232('PreLoadRaw',TDT.din,s232('DSPid',2),MONO,STACK,[],[],1,1,1);
s232('pushf',HRTF.AdaptR,PDR.HRTF_nlines);
s232('PreLoadRaw',TDT.din,s232('DSPid',3),MONO,STACK,[],[],1,1,1);
s232('dropall');

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

if(size(HRTF.TestL,1)~=1 && size(HRTF.TestL,2)~=PDR.HRTF_nlines)
    h = warndlg('Test HRTF arrays must be 1 x HRTF_nlines');
    uiwait(h);
    out=-1;
end
if(size(HRTF.TestR,1)~=1 && size(HRTF.TestR,2)~=PDR.HRTF_nlines)
    h = warndlg('Test HRTF arrays must be 1 x HRTF_nlines');
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