function speakerTestSetup_lds()
global PDR
% PARAMS

PDR.filename = 'speaker_test';
cd(PDR.data_path);
if exist([PDR.filename '.mat'])
    delete([PDR.filename '.mat']);
end
PDR.comments = 'Testing Equipment!';
PDR.record = 1;
PDR.isi_buf = 5;
% ORDERED LIST OF SPEAKER LOCATIONS
PDR.npretrials = 1;
PDR.TEST_trial_freq = 1;
PDR.n_test_trials = 2*length(PDR.TEST_azimuths);
ntrials = PDR.npretrials + PDR.TEST_trial_freq*PDR.n_test_trials;
PDR.SOUNDS_location_sequence = NaN*ones(2,ntrials+1);
PDR.SOUNDS_location_sequence(1,:) = PDR.LEAD_pos*ones(1,ntrials + 1); % lead positions in first row
PDR.SOUNDS_location_sequence(2,:) = PDR.LAG_hab_pos*ones(1,ntrials + 1); % lag positions in second row

PDR.TEST_trial_jitter=0;
TEST_trial_jitter = PDR.TEST_trial_jitter;
no_angles = length(PDR.TEST_azimuths);

% create a trial sequence and convert to location sequence 
% (2 rows, 1st row= El, 2nd row = Az)
last_test=NaN;

% create list of locations with randomized test locs
cnt = PDR.npretrials; % skip pretrials
id = 0;
while cnt < ntrials
    cnt = cnt + PDR.TEST_trial_freq;
    id = id + 1;
    PDR.SOUNDS_location_sequence(2,cnt)=PDR.TEST_azimuths(mod(id,length(PDR.TEST_azimuths))+1);
    cnt = cnt + round(TEST_trial_jitter*rand); % jitter test trial occurances (if requested)
    if cnt >= (ntrials-PDR.TEST_trial_freq)
        break;
    end
end

% setup rove sequence:
num_snds=PDR.SOUNDS_num_carriers;
PDR.SOUNDS_rove_sequence = ceil(num_snds-num_snds*rand(1,ntrials+1));
calcSessionLen_lds;
sec=num2str(PDR.len_session(2));
if length(sec)==1
    sec=['0' sec];
end
minutes=num2str(PDR.len_session(1));
if length(minutes)==1
    minutes=['0' minutes];
end
% TELL USER HOW LONG EQUIPMENT TESTING WILL TAKE!
h=msgbox(['Testing equipment will take: ' minutes ' minutes & ' sec ' seconds']);
uiwait(h);
% MAKE SOUND TOKENS:
AMStim; % calculate lead lag sounds using correlation and state values (for reproducible sounds)
soundBufferSetup_lds;
% SETUP CALIBRATED SCALES & ATTENS:
readCalibFiles_lds;