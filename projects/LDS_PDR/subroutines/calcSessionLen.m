function calcSessionLen()
global PDR
% calculate length of session:
tmp = (PDR.npretrials + PDR.TEST_trial_freq*PDR.n_test_trials)...
    *(PDR.isi_buf+1)*(PDR.buf_pts/PDR.stim_Fs); % length in seconds
PDR.len_session(1)=floor(tmp/60); % Calculate length of session (in minutes)!
PDR.len_session(2)=round(tmp-floor(tmp/60)*60); % seconds