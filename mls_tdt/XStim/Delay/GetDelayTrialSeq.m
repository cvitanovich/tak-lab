function [TrialSeq, ColumnLabels, DelayModLabels] = GetDelayTrialSeq()
    % GetDelayTrialSeq
    global H
	global XStimParams
	global TDT
	global FN
	global C_
	global M_
	global GUI
    
    randSeq = 1; % use random sequence or linear order
    
    %Get all HRTF spectra indices
	XStimParams.locations = GUI.locations1';
	if ~size(XStimParams.locations,2) % | size(XStimParams.locations,2)>1 % Limit trials to one location
        set(H.pickerfig,'Color', [.1 .2 .8]);
        set(H.picker_error,'visible','on');
        pause;
        XStimParams.locations = GUI.locations1';
        set(H.picker_error,'visible','off');
        set(H.pickerfig,'Color', [.8 .8 .8]);
	end
       
    % 50 random noise states (corresponding with saved files).
    RandStates = [...
    73 167 617 364 662 593 538 194 853 610 294 ...
    479 71 105 162 770 143 116 252 101 377 706 ...
    273 574 915 661 935 355 530 540 220 232 886 ...
    70 65 571 35 339 87 281 795 283 974 248 ...
    995 936 769 943 127 224];
    StateLow = 1; StateHi = 50; % choose state between 1 and 50
    
    Delays = XStimParams.DelayTimes(~isnan(XStimParams.DelayTimes));
    OnRamps = XStimParams.ramp_timeS_on(~isnan(XStimParams.ramp_timeS_on));
    OffRamps = XStimParams.ramp_timeS_off(~isnan(XStimParams.ramp_timeS_off));
    LagExts = XStimParams.lag_seg_ext(~isnan(XStimParams.lag_seg_ext));
    
    % Uses DelayMods to alter delay type
    % 1 =normal
    % 2 =uncorr
    % 3 =corrrel & no lead segment
    % 4 =corrrel & no lag segment
    % 5 =uncorrrel & no lead segment
    % 6 =uncorrrel & no lag segment
     DelayModLabels = {'Normal' 'Uncorrel' 'Correl no lead seg' 'Correl no lag seg' 'Uncorrel no lead seg' 'Uncorrel no lag seg' }; 
     DelayMods = [];

    % getting these here because they don't seem to be updated properly...
    XStimParams.uncorrel = get(H.uncorrel,'Value');
    uncorrel=XStimParams.uncorrel;
    XStimParams.SilentDelay = get(H.silentdelay,'Value');
    SilentOnset=XStimParams.SilentDelay;
    XStimParams.SilentLag = get(H.silentlag,'Value');
    SilentLag=XStimParams.SilentLag;
    
    if(uncorrel & SilentOnset & SilentLag)
        DelayMods = [1 2 3 4]; 
        % 1=normal, 2=uncorr, 3=corr no lead segment, 4=corr no lag segment
    else
        if(uncorrel & SilentOnset & ~SilentLag)
            %DelayMods = [1 2 3];
            DelayMods = [1 2 3 5]; 
            % 1=normal, 2=uncorr, 3=corr no lead segment, , 5=uncorr & no lead segment
        else
            if(~uncorrel & SilentOnset & ~SilentLag)
                DelayMods = [1 3]; % 3=corr without onset
            else
                if(uncorrel & ~SilentOnset & ~SilentLag)
                    DelayMods = [1 2]; % 1=normal, 2=uncorr
                else
                    if(~uncorrel & ~SilentOnset & SilentLag)
                        DelayMods = [1 4]; % 1=normal, 4=corr no lag segment
                    else
                       if(~uncorrel & ~SilentOnset & ~SilentLag)
                            DelayMods = [1]; % 1=normal
                       else
                           warning('no valid stimulus combination selected');
                           DelayMods = [1]; % 1=normal
                       end
                    end
                end
            end
        end
    end

     % NEW  METHOD...
%     if PStimG.DelayMods(1,1)
%         DelayMods = [DelayMods 1];
%     end
%     if PStimG.DelayMods(1,2)
%         DelayMods = [DelayMods 2];
%     end
%     if PStimG.DelayMods(1,3)
%         DelayMods = [DelayMods 3];
%     end
%     if PStimG.DelayMods(1,4)
%         DelayMods = [DelayMods 4];
%     end
%     if PStimG.DelayMods(1,5)
%         DelayMods = [DelayMods 5];
%     end
%     if PStimG.DelayMods(1,6)
%         DelayMods = [DelayMods 6];
%     end
%     if isempty(DelayMods)
%         DelayMods = [DelayMods 1];
%         disp('Added normal condition since none were selected...');
%     end
        
     
    % number of trials = product of each pertubation
    numtrials = length(Delays) * length(OnRamps) * length(OffRamps) * length(LagExts) * length(DelayMods);
    
%     disp(length(Delays));
%     disp(length(OnRamps));
%     disp(length(OffRamps));
%     disp(length(LagExts));
%     disp(length(DelayMods));
        
    ColumnLabels = {'Trial' 'Rep' 'Az' 'El' 'AzOffset' 'ElOffset' 'State' 'Atten' };    
    ColumnLabels = [ColumnLabels {'Delay' 'OnRamp' 'OffRamp' 'LagExt' 'DelayMod'}];
    
    numVariables = 20; % make large so that new variables can be added
    
    TrialSeq = zeros((numtrials*XStimParams.numreps) + XStimParams.numreps, numVariables);  % add PXStimParams..numreps for NaN blanks
    
    % Get ordered sequence
    DelaySeq=zeros(numtrials,1);
    OnRampSeq=zeros(numtrials,1);
    OffRampSeq=zeros(numtrials,1);
    LagExtSeq=zeros(numtrials,1);
    DelModSeq=zeros(numtrials,1);
    N = 1;
    for DelN = 1:length(Delays)
        for OnRmpN = 1:length(OnRamps)
            for OffRmpN = 1:length(OffRamps)
                for LgExtN = 1:length(LagExts)
                    for DelModN = 1:length(DelayMods)
                        DelaySeq(N)=Delays(DelN);
                        OnRampSeq(N)=OnRamps(OnRmpN);
                        OffRampSeq(N)=OffRamps(OffRmpN);
                        LagExtSeq(N)=LagExts(LgExtN);
                        DelModSeq(N)=DelayMods(DelModN);
                        N=N+1;
                    end
                end
            end
        end
    end
   
    cnt=1;
    for rep = 1:XStimParams.numreps
        % insert NaNs to trigger blank sound for baseline calculation
        TrialSeq(cnt,1) = cnt;
        TrialSeq(cnt,2) = rep;
        TrialSeq(cnt,3:size(TrialSeq,2)) = [NaN];
        cnt=cnt+1;
        
            rand('state',sum(100*clock));
            if randSeq % random sequence
                randseq = randperm(numtrials);
            else
                randseq = [1:numtrials];
            end
                      
            for perm = 1:numtrials
                TrialSeq(cnt,1) = cnt;
                TrialSeq(cnt,2) = rep;
                TrialSeq(cnt,3) = XStimParams.locations(1);
                TrialSeq(cnt,4) = XStimParams.locations(2);
                TrialSeq(cnt,5) = XStimParams.offset_az(2);
                TrialSeq(cnt,6) = XStimParams.offset_el(2);
                % rand noise state
                TrialSeq(cnt,7) = RandStates(round(StateLow+(StateHi-StateLow)*rand(1)));
                TrialSeq(cnt,8) = XStimParams.curr_ABL; 
                TrialSeq(cnt,9) = DelaySeq(randseq(perm));
                TrialSeq(cnt,10) = OnRampSeq(randseq(perm));
                TrialSeq(cnt,11) = OffRampSeq(randseq(perm));
                TrialSeq(cnt,12) = LagExtSeq(randseq(perm));
                TrialSeq(cnt,13) = DelModSeq(randseq(perm));
                cnt=cnt+1;
            end
    end
    
    %disp(size(TrialSeq,1));
    if 1
        % Eliminate trials where delay exceeds stim. duration - the duration of
        % on/offset ramps and where the onset or offset segment is silenced
        SeqLen = size(TrialSeq,1);
        %disp([1 SeqLen]);
        remove = [];
        for tr = 1:SeqLen
            if (TrialSeq(tr,13)==3 | TrialSeq(tr,13)==4 | TrialSeq(tr,13)==5 | TrialSeq(tr,13)==6) ...
                    & abs(TrialSeq(tr,9)) >= XStimParams.curr_stimdur-(max(TrialSeq(:,10))+max(TrialSeq(:,11)));
                remove=[remove tr];
            end
        end
        TrialSeq(remove,:)=[];
        if(size(TrialSeq,1)<SeqLen)
            disp(['1) Eliminated ' num2str(SeqLen-size(TrialSeq,1)) ' trials where the lead or lag segment was silenced and delay exceeded stimulus duration minus max(onset,offset) duration']);
        end

        % Eliminate trials where the lag segment is extended but is silenced
        SeqLen = size(TrialSeq,1);
        %disp([2 SeqLen]);
        remove = [];
        for tr = 1:SeqLen
            if (TrialSeq(tr,13) == 4 | TrialSeq(tr,13) == 6) & TrialSeq(tr,12) > 0;
                remove=[remove tr];
            end
        end
        TrialSeq(remove,:)=[];
        if(size(TrialSeq,1)<SeqLen)
            disp(['2) Eliminated ' num2str(SeqLen-size(TrialSeq,1)) ' trials where the lag segment was extended but silenced']);
        end

         % Eliminate trials where DelayTime < -DUR but the lag is extended or silenced
         % This might now be accomplished above???????????

        % Eliminate trials where the lag segment is extended and  delay > 10 | < 5 ms (why < 5?) 
        SeqLen = size(TrialSeq,1);
        %disp([3 SeqLen]);
        remove = [];
        for tr = 1:SeqLen
            if TrialSeq(tr,12) > 0 & (abs(TrialSeq(tr,9)) > 10 | abs(TrialSeq(tr,9)) < 5);
                remove=[remove tr];
            end
        end
        TrialSeq(remove,:)=[];
        if(size(TrialSeq,1)<SeqLen)
            disp(['3) Eliminated ' num2str(SeqLen-size(TrialSeq,1)) ' trials where the lag segment was extended but delay > 10 | < 5 ms']);
        end

        % Eliminate trials where the lag extenstion is negative but the lag
        % segment is >= 15 ms
        SeqLen = size(TrialSeq,1);
        %disp([4 SeqLen]);
        remove = [];
        for tr = 1:SeqLen
            if abs(TrialSeq(tr,9)) + TrialSeq(tr,12) >= 15 & TrialSeq(tr,12) < 0;
                remove=[remove tr];
            end
        end
        TrialSeq(remove,:)=[];
        if(size(TrialSeq,1)<SeqLen)
            disp(['4) Eliminated ' num2str(SeqLen-size(TrialSeq,1)) ' trials where the lag extention was negative but the segment was >= 15 ms']);
        end

        % Eliminate trials where NEGATIVE lag Extensions are longer than or
        % equal to lag delay
        SeqLen = size(TrialSeq,1);
        %disp([5 SeqLen]);
        remove = [];
        for tr = 1:SeqLen
            if abs(TrialSeq(tr,9)) <= -TrialSeq(tr,12) & TrialSeq(tr,12) < 0;
                remove=[remove tr];
            end
        end
        TrialSeq(remove,:)=[];
        if(size(TrialSeq,1)<SeqLen)
            disp(['5) Eliminated ' num2str(SeqLen-size(TrialSeq,1)) ' trials where a negative lag extention was longer than or equal to delay']);
        end
    end
    