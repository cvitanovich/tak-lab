function [quit,PDR] = setupTrialSeq_spatialPDR(PDR)
% function to set up trial sequence for spatialPDR


setup=1;
while setup
    % Setup randomized scale and loc sequences:
    scales = PDR.TEST_scales;
    outlier_scales = PDR.TEST_outlier_scales;
    
    mu = mean(scales);
    scales = sort(scales);
    
    % randomized test scales:
    scale_seq = zeros(1,PDR.ntrials+1);
    start = PDR.npretrials+1;
    stop = PDR.ntrials+1 - mod(PDR.ntrials-PDR.npretrials+1,PDR.TEST_trial_freq);
    nm = floor((stop - start)/PDR.TEST_trial_freq)+1;
    tmp=ceil(length(scales)*rand(1,nm));
    scale_seq(start:PDR.TEST_trial_freq:stop) = scales(tmp);

    % add outliers
    test_trial_ids=find(scale_seq~=0);
    len=length(test_trial_ids);
    for(j=1:PDR.TEST_n_outliers)
        rand_scale=outlier_scales(ceil(length(outlier_scales)*rand)); % pick a random outlier scale
        rand_id=ceil(len*rand); % pick a random test trial
        scale_seq(test_trial_ids(rand_id))= rand_scale; % replace this test trial with an outlier
    end
    
    PDR.TEST_scale_sequence = scale_seq;
    
    % create a sequence of locations that is randomized:
    start = PDR.npretrials+1;
    stop = PDR.ntrials+1 - mod(PDR.ntrials-PDR.npretrials+1,PDR.TEST_trial_freq);
    nm = floor((stop - start)/PDR.TEST_trial_freq)+1;
       
    % run script to plot trials
    screen_size = get(0, 'ScreenSize'); hFig = figure;
    set(hFig, 'Position', [0.05*screen_size(3) 0.05*screen_size(4) 0.7*screen_size(3) 0.8*screen_size(4)] );
    hold on;
    PDR=plot_trials_spatialPDR(PDR);
    
    warning off MATLAB:QUESTDLG:stringMismatch;
    test=questdlg(['# Tests = ' num2str(PDR.n_test_trials) ' ... Acceptable trial sequence?'],'Trial Sequence Confirmation','YES','NO','QUIT','NO');
    
    if strcmp(test,'YES')
        setup = 0;
    end
    if strcmp(test,'QUIT')
        quit = 1;
        close;
        return;
    end
    if strcmp(test,'NO')
        close;
    end
end
quit = 0;