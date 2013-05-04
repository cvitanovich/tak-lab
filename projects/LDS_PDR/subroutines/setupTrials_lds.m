function setupTrials_lds()
global PDR

%disp('Setting up Trial Sequence')
% SETUP randomized trial IDs
ntrials = PDR.npretrials + PDR.TEST_trial_freq*PDR.n_test_trials;
PDR.SOUNDS_location_sequence = NaN*ones(2,ntrials+1);
PDR.SOUNDS_location_sequence(1,:) = PDR.LEAD_pos*ones(1,ntrials + 1); % lead positions in first row
PDR.SOUNDS_location_sequence(2,:) = PDR.LAG_hab_pos*ones(1,ntrials + 1); % lag positions in second row

TEST_trial_jitter = PDR.TEST_trial_jitter;
no_angles = length(PDR.TEST_azimuths);

% create a trial sequence and convert to location sequence 
% (2 rows, 1st row= El, 2nd row = Az)
%PDR.n_outlier_trials = round(PDR.outlier_fraction*PDR.n_test_trials);
separations=PDR.TEST_azimuths - PDR.LAG_hab_pos;
outlier_IDs=find(abs(separations)>=PDR.outlier_separation);
not_outlier_IDs=find(abs(separations)<PDR.outlier_separation);
rep_freq=5; % don't repeat an outlier speaker location less than 7 test trials apart
rand_list = randomized_trials_lds(not_outlier_IDs,outlier_IDs,PDR.n_outlier_trials,PDR.n_test_trials,rep_freq);

last_test=NaN;
if length(rand_list) == 0
    return;
end
% create list of locations with randomized test locs
cnt = PDR.npretrials; % skip pretrials
id = 0;
while cnt < ntrials
    cnt = cnt + PDR.TEST_trial_freq;
    id = id + 1;
    PDR.SOUNDS_location_sequence(2,cnt)=PDR.TEST_azimuths(rand_list(id));
    cnt = cnt + round(TEST_trial_jitter*rand); % jitter test trial occurances (if requested)
    if cnt >= (ntrials-PDR.TEST_trial_freq)
        break;
    end
end

% setup rove sequence:
num_snds=PDR.SOUNDS_num_carriers;
PDR.SOUNDS_rove_sequence = ceil(num_snds-num_snds*rand(1,ntrials+1));
