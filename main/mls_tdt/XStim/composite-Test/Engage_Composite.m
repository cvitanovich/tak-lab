function [] = Engage_Composite(FNlist, Test_array)

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI
global F

%Engage_Composite
% Initialize Application and get AP2 and XBUS locks
if(S232('S2init', 0, 'INIT_PRIMARY', 20000) == 0)
   disp('Cannot initialize a primary process')
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

% plays out stimulus pairs saved to temp_stim_path
stimuli_dir = FN.temp_stim_path;
fclose all;

XStimParams.curr_stimdur = str2num(get(H.composite_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.composite_ISI,'String'));
XStimParams.numreps = str2num(get(H.composite_numreps,'String'));
XStimParams.reset_flag = 0;

%Specify DAMA buffers
clear BUF
BUF.L1				= 1;
BUF.R1				= 2;
BUF.playseq_L1		= 5;
BUF.playseq_R1		= 6;
BUF.playspec1		= 7;

%Make play sequence buffers
S232('allot16',BUF.playseq_L1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.L1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_L1);

S232('allot16',BUF.playseq_R1,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.R1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_R1);

%Make play specification buffer
S232('allot16',BUF.playspec1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L1);
S232('make',1,BUF.playseq_R1);
S232('make',2,0);
S232('qpop16',BUF.playspec1);

% Add a piece of silence prior to stimulus to calculate spontaneous rate
DUR_silence = XStimParams.silence_lead; 					%ms
silence_len = (DUR_silence * round(TDT.Fs/1000));
% Add a piece of silence after stimulus 
DUR_silence2 = XStimParams.silence_trail; 					%ms
silence_len2 = (DUR_silence2 * round(TDT.Fs/1000));

%Make Stimulus buffers
DUR = XStimParams.curr_stimdur;
Npts_totalplay = (DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen;

S232('allot16',BUF.L1,Npts_totalplay);
S232('allot16',BUF.R1,Npts_totalplay);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,Npts_totalplay);

%Get Earphone filters
if FN.HRTFfiletype(6) == 1
    ephonefname = [FN.ephone_path FN.ephone2];
    ephonefilt_left  = (mtlrch(ephonefname,1))';
    ephonefilt_right = (mtlrch(ephonefname,2))';
else
    eval(['load -mat ' FN.ephone_path FN.ephone2]);
    ephonefilt_left  = TF1;
    ephonefilt_right = TF2;
    clear TF1 TF2 dir
end

%Load Earphone filters
dspid_left = 0; dspid_right = 1;
S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('FFF'));           % used to be '0xFFF'
S232('dropall');
%Make connections for left ear
S232('PD1addsimp',1,S232('DSPout',dspid_left),S232('DAC',0)); %DSPout to DAC0
S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left)); %IB to DSPin
%Make connections for right ear
S232('PD1addsimp',1,S232('DSPout',dspid_right),S232('DAC',1));
S232('PD1specIB',1,S232('IB',1),S232('DSPin',dspid_right));
%Load left      
S232('pushf',ephonefilt_left,length(ephonefilt_left));
S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.ephonescale,1.0,1);
%Load right
S232('pushf',ephonefilt_right,length(ephonefilt_right));
S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.ephonescale,1.0,1);

set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

ISI = XStimParams.test_ISI;
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000)); 				%correct for ITD filtlength
ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));			%correct for ephonefilt

% update display
setInfo_Composite;

%Begin playing sounds   
set(H.composite_status,'String','Status: playing stimuli');
set(H.composite_status,'BackgroundColor','green');
set(H.composite_remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
buffcycle = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.composite_fig'))
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
numstims = length(FNlist);
finalspikematrix = zeros(length(Test_array),numstims);
% to convert playlist.type to Test_array.type
ind2Testarray = ones(12,1)*NaN;
for i = 1:length(Test_array)
    ind2Testarray(Test_array(i).type) = i;
end
% create list(s) of itds and ilds to be played
if ~isnan(ind2Testarray(1))
    temp_params = Test_array(ind2Testarray(1)).params;
    %inc = (temp_params.hiabl - temp_params.loabl)/(temp_params.numabls-1);
    %abls = temp_params.loabl : inc : temp_params.hiabl;
    abls = temp_params.abls;
end
if ~isnan(ind2Testarray(2))
    temp_params = Test_array(ind2Testarray(2)).params;
    %inc = (temp_params.hiitd - temp_params.loitd)/(temp_params.numitds-1);
    %itds = temp_params.loitd : inc : temp_params.hiitd;
    itds = temp_params.itds;
end
if ~isnan(ind2Testarray(3))
    temp_params = Test_array(ind2Testarray(3)).params;
    %inc = (temp_params.hiild - temp_params.loild)/(temp_params.numilds-1);
    %ilds = temp_params.loild : inc : temp_params.hiild;
    ilds = temp_params.ilds;
end

%loop for reps
while (exist1('H.composite_fig') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli:
    %       stimuli come in alphabetized (as FNlist), then are re-randomized (into playlist)
    %       the ilds, itds, abls vectors are in min-to-max order
    %       temp1 (below) provides an index into ilds, itds or abls derived
    %           from the filename being played out
    %       param1 (and param2) are created above as a matrix of NaNs
    %       below, param1 is filled with the indices into ilds, itds or
    %           abls (from temp1), OR param1 and param2 are filled with the
    %           EL and AZ, respectively, from the filename being played out
    %       the ordering of param1 and param2 is that of FNlist, since
    %       playlist is re-randomized for each rep

    % initialize param1,2 for data storage
    param1 = ones(length(Test_array),length(FNlist))*NaN;
    param2 = ones(length(Test_array),length(FNlist))*NaN;

    randseq = randperm(numstims);   
    
    playlist = FNlist(randseq);
    tic
    for trialnum = 1:numstims
        %Check for pause by user
        if pause_check(length(Test_array))    return; end
        clear temp*
        ind_type = ind2Testarray(playlist(trialnum).type);         % index into Test_array for this stim
        % get el, az, ild, or itd values and place into param1 and param2
        ind2name = findstr(playlist(trialnum).name,'.');
        ind2name = findstr(playlist(trialnum).name(ind2name:end),'_')+ind2name;
        switch length(ind2name)
            case 1
                if findstr(playlist(trialnum).name,'ITD')
                    temp1 = find(itds ==str2num(playlist(trialnum).name(ind2name(1):end)));  % gives index to itd
                elseif findstr(playlist(trialnum).name,'ILD')
                    temp1 = find(ilds ==str2num(playlist(trialnum).name(ind2name(1):end)));  % gives index to ild
                elseif findstr(playlist(trialnum).name,'ABL')
                    temp1 = find(abls ==str2num(playlist(trialnum).name(ind2name(1):end)));  % gives index to abl
                elseif findstr(playlist(trialnum).name,'itd')
                    temp1 = find(itds ==str2num(playlist(trialnum).name(ind2name(1):end)));  % gives index to itd
                elseif findstr(playlist(trialnum).name,'ild')
                    temp1 = find(ilds ==str2num(playlist(trialnum).name(ind2name(1):end)));  % gives index to ild
                elseif findstr(playlist(trialnum).name,'abl')
                    temp1 = find(abls ==str2num(playlist(trialnum).name(ind2name(1):end)));  % gives index to abl
                end
            case 2
                temp1 = str2num(playlist(trialnum).name(ind2name(1):ind2name(2)-2));
                temp2 = str2num(playlist(trialnum).name(ind2name(2):end));
                param2(ind_type,randseq(trialnum)) = temp2;
        end
        param1(ind_type,randseq(trialnum)) = temp1;
        
         if playlist(trialnum).type ==1
            S232('PA4atten',1,abs(abls(temp1))-20);					% HB1 adds 20 dB attenuation
            S232('PA4atten',2,abs(abls(temp1))-20);
        else
            S232('PA4atten',1,abs(ABL)-20);
            S232('PA4atten',2,abs(ABL)-20);
        end

        % load stimuli
        % re-zero and load left buffer
        S232('dpush',Npts_totalplay);
        S232('value',0);
        S232('qpop16',BUF.L1);
        fname = playlist(trialnum).name;
        set(H.composite_status,'String',['playing: ' fname]);
        evalstr = strcat(['S232(''pushdisk16''' ','], [' ''' FN.temp_stim_path], ...
            fname,[ ''' '], ');');
        eval(evalstr);
        evalstr = ['S232(''qpoppart16'''  ',BUF.L1,0);'];
        eval(evalstr);
        % re-zero and load right buffer
        temp = findstr('left',fname);
        fname = [fname(1:temp-1) 'right' fname(temp+4:end)];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
        eval(evalstr);
        evalstr = ['S232(''qpoppart16'''  ',BUF.R1,0);'];
        S232('dpush',Npts_totalplay);
        S232('value',0);
        S232('qpop16',BUF.R1);
        eval(evalstr);
        
        S232('seqplay',BUF.playspec1);
        S232('PD1arm',1);
        
        %Send trigger
        %Set up MII
        m100x( C_.INIT );
        m110dx( C_.INIT );
        m110dx( C_.CLOCK, mii_us_per_sample);
        m110dx( C_.MODE, M_.PST );
        
        while toc < ISI/1000 
            if pause_check(length(Test_array))    return; end
        end
        
        %Start clock
        m110dx( C_.START);
        %Send pulse: PD1 GO!
        m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
        tic
        
        %while (S232('APactive')) pause(DUR/1000+.1); end
        pause(DUR/1000+.05);     % to bypass PD1status bug
        %while(S232('PD1status',1)) usec_delay(1000);  end
        S232('PD1stop',1);
        
        %Stop the m110 and get spikes
        m110dx( C_.STOP);
        spikes = m110dx( C_.DATA, XStimParams.curr_stimdur); 			% Take XStimParams.curr_stimdur spikes max
        ind = find(spikes ~= 0); 						% Get clock events that are spikes
        spikes = spikes(ind);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
        if ~isempty(spikes) 
            spikes = spikes(ind);
            Test_array(ind_type).spikes = ...
                [Test_array(ind_type).spikes; spikes/(1000/mii_us_per_sample)]; 
            Test_array(ind_type).repnum = ...
                [Test_array(ind_type).repnum; repnum * ones(size(spikes))]; 
            Test_array(ind_type).stimnum = ...
                [Test_array(ind_type).stimnum; randseq(trialnum) * ones(size(spikes))]; 
            Test_array(ind_type).Nspikes = ...
                [Test_array(ind_type).Nspikes; length(spikes) * ones(size(spikes))]; 
            Test_array(ind_type).param1 = ...            
                [Test_array(ind_type).param1; temp1 * ones(size(spikes))]; 
            if length(ind2name)==2
                Test_array(ind_type).param2 = ...
                    [Test_array(ind_type).param2; temp2  * ones(size(spikes))];
            end
        end
        
        finalspikematrix(ind_type,randseq(trialnum)) = finalspikematrix(ind_type,randseq(trialnum)) + length(spikes);
        if pause_check(length(Test_array))    return; end
        
        remtrials = numstims - trialnum;
        set(H.composite_remtrials,'String',num2str(remtrials));
        pause(0);
    end %end loop over trials
    
    %Plot
    for i = 1:length(Test_array)
        ind = find(~isnan(param1(i,:)));            % index to param1,2 for this test
        if isempty(Test_array(i).param2)
            switch Test_array(i).type
                case 1
                    temp = abls(param1(i,ind));
                case 2
                    temp = itds(param1(i,ind));
                case 3
                    temp = ilds(param1(i,ind));
            end
            plot_stuff(Test_array(i), i, finalspikematrix(i,ind)/repnum, temp);    
        else
            plot_stuff(Test_array(i), i, finalspikematrix(i,ind)/repnum, param1(i,ind), param2(i,ind));
        end
    end
    
    %Record Data
    if(exist1('H.composite_fig') )
        for i = 1:length(Test_array)
            if isempty(Test_array(i).param2)
                datamatrix = [Test_array(i).Nspikes Test_array(i).spikes Test_array(i).repnum Test_array(i).param1];
            else                
                datamatrix = [Test_array(i).Nspikes Test_array(i).spikes Test_array(i).repnum Test_array(i).param1 Test_array(i).param2];
            end
            FN.data = Test_array(i).FN;
            temp = Test_array(i);
            record_data3(Test_array(i).params,datamatrix);
            eval(['save ' FN.data_path FN.data '.mat playlist temp -append']);
        end
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.composite_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

S232('PA4mute',1);
S232('PA4mute',2);

set(H.composite_status,'String','Status: results');
set(H.composite_status,'BackgroundColor','blue');

set(H.composite_exit,'Visible','on');
set(H.composite_reset,'Visible','on');

update_dataFN;

%%%%%%%%%
function [flag] = pause_check(Test_num)

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI
global F

flag = 0;
%Check for pause by user
while (exist1('H.composite_fig') & get(H.composite_pause,'Value'))
    pause(0);
    if(~exist1('H.composite_fig')) return; end         
    set(H.composite_exit,'Visible','on');
    set(H.composite_reset,'Visible','on');
    if(exist1('H.composite_fig') & get(H.composite_reset,'Value') == 1)
        set(H.composite_reset,'Value',0);
        set(H.composite_pause,'Value',0);
        Reset_Composite;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end

%%%%%%%%%%%
function [] = plot_stuff(Test_array, nFig, nspikes, param1, param2)

global H

switch Test_array.type
    case 1      % ABL
        [param1,ind] = sort(param1);
        nspikes = nspikes(ind);
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','ABL test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig));
        plot(param1,nspikes,'g-*','LineWidth',1.5);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('ABL (-dB atten)'); ylabel('Spike Rate (spikes/stim)');
        title('ABL Test');
    case 2      % ITD
        [param1,ind] = sort(param1);
        nspikes = nspikes(ind);
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','ITD test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig));
        plot(param1,nspikes,'g-*','LineWidth',1.5);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('ITD (us)'); ylabel('Spike Rate (spikes/stim)');
        title('ITD Test');
    case 3      % ILD
        [param1,ind] = sort(param1);
        nspikes = nspikes(ind);
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','ILD test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig));
        plot(param1,nspikes,'g-*','LineWidth',1.5);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('ILD (dB)'); ylabel('Spike Rate (spikes/stim)');
        title('ILD Test');
    case 4      % space1 FC
        dir(1,:) = param1;   
        dir(2,:) = param2;   
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','fully cued space1 test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig))
        plotdiam1(dir, nspikes);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title('space1 FC');
        colorbar
    case 5      % space1 ILA
        dir(1,:) = param1;   
        dir(2,:) = param2;   
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','ILA space1 test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig))
        plotdiam1(dir, nspikes);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title('space1 ILA');
        colorbar
    case 6      % space1 ITA
        dir(1,:) = param1;   
        dir(2,:) = param2;   
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','ITA space1 test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig))
        plotdiam1(dir, nspikes);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title('space1 ITA');
        colorbar
    case 7      % space2 FC
        dir(1,:) = param1;   
        dir(2,:) = param2;   
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','fully cued space2 test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig))
        plotdiam1(dir, nspikes);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title('space2 FC');
        colorbar
    case 8      % space2 ILA
        dir(1,:) = param1;   
        dir(2,:) = param2;   
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','ILA space2 test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig))
        plotdiam1(dir, nspikes);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title('space2 ILA');
        colorbar
    case 9      % space2 ITA
        dir(1,:) = param1;   
        dir(2,:) = param2;   
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','ITA space2 test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig))
        plotdiam1(dir, nspikes);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title('space2 ITA');
        colorbar
    case 10     % altIR FC
        dir(1,:) = param1;   
        dir(2,:) = param2;   
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','fully cued altIR test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig))
        plotdiam1(dir, nspikes);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title('altIR FC');
        colorbar
    case 11     % altIR ILA
        dir(1,:) = param1;   
        dir(2,:) = param2;   
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','ILA altIR test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig))
        plotdiam1(dir, nspikes);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title('altIR ILA');
        colorbar
    case 12     % altIR ITA
        dir(1,:) = param1;   
        dir(2,:) = param2;   
        if ~exist1('H.composite_datafig') | length(H.composite_datafig) < nFig
            H.composite_datafig(nFig) = figure('Position',[700 20 550 500],...
                'Name','ITA altIR test',...
                'NumberTitle','off');
            H.composite_spikeaxes(nFig) = axes;
        end
        figure(H.composite_datafig(nFig))
        plotdiam1(dir, nspikes);
        set(H.composite_spikeaxes(nFig),'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title('altIR ITA');
        colorbar
end