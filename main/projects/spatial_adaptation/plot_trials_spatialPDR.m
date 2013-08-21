function PDR=plot_trials_spatialPDR(PDR)
% PLOT TRIALS
tmp = find(PDR.TEST_scale_sequence~=0);
if length(tmp)<1
    warndlg('Test trial sequence init fail!!!','No Trials');
    quit = 1;
    close(hFig);
    return;
end
cnt=0;

load([PDR.CALIB_directory PDR.CALIB_fname]);
eval(['mTest=CALIB_PDR.' PDR.TEST_soundtype '_COEFS(1,1);']);
eval(['bTest=CALIB_PDR.' PDR.TEST_soundtype '_COEFS(2,1);']);
t_trials=find(PDR.TEST_scale_sequence~=0);
scales=PDR.TEST_scale_sequence(t_trials);
SPLs=mTest.*log10(scales)-PDR.base_atten+bTest;
%sz=ceil(sqrt(scales)); %ceil(log10(scales)).^3; % sizes for trial markers
loc=PDR.TEST_loc(2);
scatter(t_trials,SPLs);
PDR.n_test_trials = cnt;
xlabel('Trial #');
ylabel('SPLs');
axis([0 PDR.ntrials+1 -30 90]);
whitebg(gcf,'k');
ntests=length(find(PDR.TEST_scale_sequence~=0));
title(['Trial Sequence ( ' num2str(ntests) ' test trials )']);