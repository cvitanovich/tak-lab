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
scale_list=sort(PDR.TEST_scales);
t_trials=find(PDR.TEST_scale_sequence~=0);
scales=PDR.TEST_scale_sequence(t_trials);
sz=ceil(log10(scales)).^3; % sizes for trial markers
loc=PDR.TEST_loc(2);
scatter(t_trials,loc.*ones(1,length(t_trials)),sz,'filled')

PDR.n_test_trials = cnt;
xlabel('Trial #');
ylabel('Location (az)');
y_tick_label = ['(' num2str(PDR.TEST_loc(1)) ',' num2str(PDR.TEST_loc(2)) ')'];
set(gca,'YTick',loc,'YTickLabel',y_tick_label);
axis([0 PDR.ntrials+1 -100 100]);
whitebg(gcf,'k');
title('Test Trial Scale Sequence');