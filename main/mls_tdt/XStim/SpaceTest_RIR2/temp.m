%%%%% Loop to make the stimuli we'll use
% numtrials, trialnum reflects all iFiles and all iLocs
remreps = 1;
set(H.space_RIR_status,'String','Status: building stimuli');
set(H.space_RIR_status,'BackgroundColor','yellow');
set(H.space_RIR_status,'ForegroundColor','blue');
set(H.space_RIR_remreps,'String',num2str(remreps));
repnum = 1;
nLocs = size(XStimParams.locations,2);
nFiles = size(FN.RIR,2);
numtrials = nLocs * nFiles;
finalspikematrix = zeros(1,numtrials);

%Randomize stimFNs
FN.RIR = FN.RIR(randperm(nFiles));

% clear out temporary directory
delete( [FN.temp_stim_path '\*.*'])

trialnum = 1;
for iFile = 1:nFiles
    %Check for pause by user
    if pause_check    return; end
    % get stimulus from file
    if strcmp('mat',FN.RIR{iFile}(end-2:end))
        eval(['load ' FN.stim_path2 FN.RIR{iFile}])
        trial_left  = stim;
    else
        fid = fopen([FN.stim_path2 FN.RIR{iFile}],'r');
        trial_left = fread(fid,inf,'float');
        fclose(fid);
    end
    trial_left = trial_left(:)' - mom(trial_left,1);
    
    % trailing pad stim length of short stims
    if length(trial_left)/TDT.Fs *1000 < XStimParams.curr_stimdur
        npts = ceil(XStimParams.curr_stimdur*30 - length(trial_left));
        trial_left = [trial_left zeros(1,npts)];
    end
    
    % normalize to ACPower
    if ~XStimParams.HiDynamicRange
        trial_left = trial_left / mom(trial_left,2);
    end
    trial_right = trial_left;
    trial_left0 = trial_left;
    trial_right0 = trial_right;
    
    for iLoc = 1:nLocs        %Apply HRTF filtering
        if pause_check    return; end        %Check for pause by user
        
        if FN.HRTFfiletype(1,2) == 1
            hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locind(iLoc))-1);
            hrtf_right = mtlrch([FN.space_path FN.space_std],2*locind(iLoc));
        else
            hrtf_left = TF1_space(locind(iLoc),:);
            hrtf_right = TF2_space(locind(iLoc),:);
        end
        str1 = 'FC';
        
        
        trial_left = conv(trial_left0, hrtf_left);
        trial_right = conv(trial_right0, hrtf_right);
        
        % remove DC offset
        trial_left = trial_left - mean(trial_left);
        trial_right = trial_right - mean(trial_right);
        
        % adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
        % presented through DSPs (*.eq) at 0,0 (broadband)
        trial_left = trial_left * TDT.hrtf_Lfactor;
        trial_right = trial_right * TDT.hrtf_Rfactor;
        
        %scale stimuli 3/7/07
        % note scale factor changes with XStimParams.HiDynamicRange
        trial_left = trial_left * scaleFactor;
        trial_right = trial_right * scaleFactor;
        
        %Add in the leading silent period
        trial_left =  [zeros(1,silence_ptsLEAD) trial_left];
        trial_right = [zeros(1,silence_ptsLEAD) trial_right];
        
        %Add in the trailing silent period
        trial_left =  [trial_left zeros(1,silence_ptsLAG)];
        trial_right = [trial_right zeros(1,silence_ptsLAG)];
        
        %pad with zeros for ephonefilters
        trial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
        trial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
        
        if length(trial_left) ~= Npts_totalplay | length(trial_right) ~= Npts_totalplay
            error('engage499: filesize incorrect')
        end
        
        % save stimuli to disk
        S232('push16',trial_left,length(trial_left));
        S232('qpop16',BUF.L1);
        fname = ['RIR' num2str(iFile) '_' str1 '.left_' ...
            num2str(hrtfdirmat(1,locind(iLoc))) ...
            '_' num2str(hrtfdirmat(2,locind(iLoc)))];
        evalstr = ['S232(''dama2disk16'',BUF.L1,' [' ''' FN.temp_stim_path fname ''' '] ',0);'];
        eval(evalstr);
        
        S232('push16',trial_right,length(trial_right));
        S232('qpop16',BUF.R1);
        fname = ['RIR' num2str(iFile) '_' str1 '.right_' ...
            num2str(hrtfdirmat(1,locind(iLoc))) ...
            '_' num2str(hrtfdirmat(2,locind(iLoc)))];
        evalstr = ['S232(''dama2disk16'',BUF.R1,' [' ''' FN.temp_stim_path fname ''' '] ',0);'];
        eval(evalstr);
        
        % set up for next trial
        remtrials = numtrials - trialnum;
        set(H.space_RIR_remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end     % iLoc
end 	% iFile

%%%%%%%%%%%%%%%%%%% finished making stimuli