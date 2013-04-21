%*******************************************************************************
%	The SEARCH operation
%*******************************************************************************

% Initialize Application and get AP2 and XBUS locks
if(S232('S2init', 0, 'INIT_PRIMARY', 10000) == 0)
    disp('Cannot initialize a secondary process')
    return;
end

if(S232('APlock', 100, 0) == 0)
    disp('Cannot acquire lock on AP2 Card')
    s232('S2close');
    return;
end

if(S232('XBlock', 100, 0) == 0)
    disp('Cannot acquire lock on X Bus')
    s232('APunlock', 0);
    s232('S2close');
    return;
end

%Set up DAMA buffers
%Specify DAMA buffers
BUF.left_1			= 1;
BUF.right_1		= 2;
BUF.playseq_left_1			= 3;
BUF.playseq_right_1		= 4;
BUF.playspecbuf_1			= 5;
BUF.left_2			= 6;
BUF.right_2		= 7;
BUF.playseq_left_2			= 8;
BUF.playseq_right_2		= 9;
BUF.playspecbuf_2			= 10;

testbuf_left			= 11;
testbuf_right			= 12;

%Make Stimulus buffers
DUR = get(H.DUR,'Value');
S232('allot16',BUF.left_1,(DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
S232('allot16',BUF.right_1,(DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
S232('allot16',BUF.left_2,(DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
S232('allot16',BUF.right_2,(DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
S232('allot16',testbuf_left,(DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
S232('allot16',testbuf_right,(DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);

%Make play sequence buffers
S232('allot16',BUF.playseq_left_1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.left_1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_left_1);

S232('allot16',BUF.playseq_right_1,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.right_1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_right_1);

S232('allot16',BUF.playseq_left_2,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.left_2);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_left_2);

S232('allot16',BUF.playseq_right_2,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.right_2);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_right_2);


%Make play specification buffer
S232('allot16',BUF.playspecbuf_1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_left_1);
S232('make',1,BUF.playseq_right_1);
S232('make',2,0);
S232('qpop16',BUF.playspecbuf_1);

S232('allot16',BUF.playspecbuf_2,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_left_2);
S232('make',1,BUF.playseq_right_2);
S232('make',2,0);
S232('qpop16',BUF.playspecbuf_2);

%Check to see if earphone filters are being used
ephoneflag = 0;
if(exist1('H.fig') & get(H.ephoneuseit,'Value'))
    if isempty(FN.ephone)
        EphoneFileDiagBox;
    end
    ephoneflag = 1;
    ephonefilt_left  = (mtlrch([FN.ephone_path FN.ephone],1))';
    ephonefilt_right = (mtlrch([FN.ephone_path FN.ephone],2))';
end

set(H.ephonefile,'Enable','off');
set(H.ephoneuseit,'Enable','off');

%Set MII parameters
mii_us_per_sample = 10; %microsecond per sample
mii_separation = 100; %only take events separated by 100 samples (i.e., 1 ms)

% Cycle over stimuli
cycle = 0;
while(exist1('H.searchfig'))
    while (exist1('H.searchfig') & get(H.pausesearch,'Value') == 1)
        pause(0);
        if(~exist1('H.searchfig')) return; end         
    end
    if(exist1('H.searchfig') & get(H.pausesearch,'Value') == 0)
        if(~exist1('H.searchfig')) return; end         
    end
    
    % stimulus-type toggle
    if get(H.togglestim,'value')
        set(H.togglestim,'string','BBnoise');
        set(H.stim_type,'Value',1);
    else
        set(H.togglestim,'string','tone');
        set(H.stim_type,'Value',8);
    end        
    
    
    if(get(H.stim_type,'Value') == 1)
        set(H.tonal_frequency_title,'Visible','on');
        set(H.tonal_frequency,'Visible','on');
        set(H.tonal_frequency_text,'Visible','on');
        set(H.tonal_frequency_hz,'Visible','on');
    else
        set(H.tonal_frequency_title,'Visible','off');
        set(H.tonal_frequency,'Visible','off');
        set(H.tonal_frequency_text,'Visible','off');
        set(H.tonal_frequency_hz,'Visible','off');
    end
    
    if(get(H.plotpsd,'Value'))
        set(H.plotpsd,'BackgroundColor','green');
        plotpsdflag = 1;
    else
        set(H.plotpsd,'BackgroundColor','blue');
        plotpsdflag = 0;
    end
    
    ILD = round(get(H.ILDslider,'Value'));
    ISI = XStimParams.search_ISI;
    ABL = round(get(H.ABLslider,'Value'));
    set(H.ILDslider_text,'String',[num2str(round(ILD)) ' dB']);
    set(H.ILDslider,'Value',(round(ILD)));
    set(H.ABLslider_text,'String',[num2str(round(ABL)) ' dB']);
    set(H.ABLslider,'Value',(round(ABL)));
    
    if(DUR ~= XStimParams.curr_stimdur) | mod(cycle,10) == 0;
        DUR = XStimParams.curr_stimdur;
        S232('deallot',BUF.left_1);
        S232('deallot',BUF.right_1);
        S232('deallot',testbuf_left);
        S232('deallot',testbuf_right);
        
        S232('allot16',BUF.left_1,(DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
        S232('allot16',BUF.right_1,(DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
        S232('allot16',testbuf_left,(DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
        S232('allot16',testbuf_right,(DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
    end
    
    %Create Stimulus
    switch get(H.stim_type,'Value')
        case 1, %Tones
            %Make tone at specified frequency and ild - NOT GETTING ILD by ATTENUATORS
            %Lright - Lleft = 20*log10(Ampright/Ampleft), Hartmann p. 29
            FREQ = get(H.tonal_frequency,'Value');
            stim = MakeTone(TDT.Fs,FREQ,DUR);
        case 8, %BroadBand Noise
            stim = MakeBBNoise(TDT.Fs,DUR);
        otherwise,
            disp('This choice of stimulus type is not supported in search mode')
            set(H.stim_type,'Value',8);
    end
    
    % apply ILD half to each side
    trial_left  = stim * (10^(-ILD/(2*20)));
    trial_right = stim * (10^(ILD/(2*20)));
    
    % 3/7/07 multiply by freefield_Xfactors
    trial_left = trial_left * TDT.freefield_Lfactor;
    trial_right = trial_right * TDT.freefield_Rfactor;
    
    % 4/10/07 multiply by .6 to make more similar to space... not sure why
    % this is, but ...
    trial_left = trial_left * .6;
    trial_right = trial_right * .6;

    
    
    % modulate stimulus
    if ~strcmp(XStimParams.mod_type,'None')
        Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
        trial_right = trial_right .* Envelope(:)';
        trial_left = trial_left .* Envelope(:)';
    end
    
    %Ramp
    ramp_time = 5; %ms
    while (S232('APactive')) pause(0); end
    [trial_left] = ramp_sound(trial_left,TDT.Fs,ramp_time);
    [trial_right] = ramp_sound(trial_right,TDT.Fs,ramp_time);
    
    %Apply ITD filtering
    itdleft = 0; itdright = 0;
    ITD = round(get(H.ITDslider,'Value'));
    if(ITD < 0)
        itdleft = 0;
        itdright = abs(ITD);
    elseif(ITD > 0)
        itdleft = abs(ITD);
        itdright = 0;
    end
    eval(['load ' FN.itd_path 'itdfilt' num2str(itdleft)]);
    eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
    clear itd_filt;
    eval(['load ' FN.itd_path 'itdfilt' num2str(itdright)]);
    eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
    clear itd_filt;
    trial_left = conv(trial_left,itd_filt_left);
    trial_right = conv(trial_right,itd_filt_right);
    
% remove DCoffset
    trial_left = trial_left - round(mean(trial_left));
    trial_right = trial_right - round(mean(trial_right));
   
    %%%%%   scaling 3/7/07
    trial_left = trial_left * TDT.scaleFactor;
    trial_right = trial_right * TDT.scaleFactor;
    
    %Wait till PD1 is finished
    %while S232('PD1status',1) usec_delay(1000); end
    usec_delay((DUR+100)*500);

    S232('PD1stop',1);
    
    %Load stimulus buffers & build play sequence and specification buffers
    S232('push16',trial_left,length(trial_left));
    S232('qpop16',BUF.left_1);
    
    S232('push16',trial_right,length(trial_right));
    S232('qpop16',BUF.right_1);
    
    %Set Attenuators
    latten = abs(ABL)-20;
    ratten = abs(ABL)-20;
    S232('PA4atten',1,latten); S232('PA4atten',2,ratten);
    
    %Stop the m110 and get spikes
    pause(0.1)									% get some spikes after stimulus
    m110dx( C_.STOP);
    spikes = m110dx( C_.DATA, 200);		 %Take 200 spikes max
    ind = find(spikes ~= 0); 				%Get clock events that are spikes
    spikes = spikes(ind);
    ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
    if(~isempty(spikes)) spikes = spikes(ind); end   
    
    %Plot the spikes
    if(exist1('H.searchfig') & ~exist1('H.spikefig'));
        H.spikefig = figure('Units','pixels',...
            'Position',[950 700 300 250],...
            'Name','RASTER PLOT');
        H.spikeaxes = axes('Position',[0.2 0.2 .7 .7]);
        plotcycle = 0;
        hold on
    end
    if(exist1('H.searchfig') & exist1('H.spikeaxes'))
        axes(H.spikeaxes);
        currXlim = get(gca,'Xlim');
        if get(H.clearplot_toggle,'value')
            cycle = 0;
            set(H.clearplot_toggle,'value',0)
        end
        plotcycle = mod(cycle,20);
        if(plotcycle == 0 | (round(length(trial_left)/(TDT.Fs/1000)+150) ~= currXlim(2)))
            set(H.spikeaxes,...
                'YDir','reverse',...
                'Color','black',...
                'Ylim',[1 20],...
                'Xlim',[0 round(length(trial_left)/(TDT.Fs/1000) + 150)],...
                'XLabel',text('String','Time (ms)'));
            T = get(gca,'children');
            if size(T,1)>1
                delete(T([1:end-1]))
            end
            plot([round(length(trial_left)/(TDT.Fs/1000)) round(length(trial_left)/(TDT.Fs/1000))],...
                [1 100],'r','LineWidth',2);
        end
        plot(spikes/(1000/mii_us_per_sample),(plotcycle+1)*ones(length(spikes)), 'c.', 'MarkerSize',5);
    end
    
    if(exist1('H.searchfig'))
        
        S232('PD1clear',1);
        S232('PD1srate', 1, 1e6/TDT.Fs);
        S232('PD1npts',1,(DUR*(round(TDT.Fs/1000))) + TDT.itdfiltlen);
        
        %Load Earphone filters
        if(get(H.ephoneuseit,'Value')) %Use earphone filters
            dspid_left = 0; dspid_right = 1;
            S232('PD1clrsched',1);
            S232('PD1nstrms',1,2,0);
            S232('PD1resetDSP',1,hex2dec('FFF'));
            S232('dropall');
            %Make connections for left ear
            S232('PD1addsimp',1,S232('DSPout',dspid_left),S232('DAC',0)); %DSPout to DAC0
            S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left)); %IB to DSPin
            %Make connections for right ear
            S232('PD1addsimp',1,S232('DSPout',dspid_right),S232('DAC',1));
            S232('PD1specIB',1,S232('IB',1),S232('DSPin',dspid_right));
            %Load left      
            while (exist1('H.searchfig') & S232('APactive'))  pause(0); end
            S232('pushf',ephonefilt_left,length(ephonefilt_left));
            while (exist1('H.searchfig') & S232('APactive'))  pause(0); end
            S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.ephonescale,1.0,1);
            %Load right
            S232('pushf',ephonefilt_right,length(ephonefilt_right));
            while (exist1('H.searchfig') & S232('APactive'))  pause(0); end
            S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.ephonescale,1.0,1);
            while (exist1('H.searchfig') & S232('APactive'))  pause(0); end
        else
            S232('PD1clrsched',1);
            S232('PD1nstrms',1,2,0);
            S232('PD1specIB',1,S232('IB',0),S232('DAC',0));
            S232('PD1specIB',1,S232('IB',1),S232('DAC',1));
        end
        
        S232('seqplay',BUF.playspecbuf_1);
        S232('PD1arm',1);
        
        %Send trigger
        %Set up MII
        m100x( C_.INIT );
        m110dx( C_.INIT );
        m110dx( C_.CLOCK, mii_us_per_sample);
        m110dx( C_.MODE, M_.PST );
        
        pause(ISI/1000);          
        %Start clock
        m110dx( C_.START);
        %Send pulse: PD1 GO!
        m101x( C_.DATA,M_.BIT,M_.PULSE,0); 				%Use port 0 for the pulse
    end
        
    cycle = cycle + 1;
    pause(0);
end