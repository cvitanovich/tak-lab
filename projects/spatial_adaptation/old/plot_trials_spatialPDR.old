function PDR=plot_trials_spatialPDR(PDR)
% PLOT TRIALS
x=0;
tmp = find(PDR.TEST_loc_sequence~=0);
if length(tmp)<1
    warndlg('No trials with this config!','No Trials');
    quit = 1;
    close(hFig);
    return;
end
cnt=0;
scale_list=sort(PDR.TEST_scales);
for i=1:length(PDR.TEST_loc_sequence)
    x=x+1;
    if PDR.TEST_loc_sequence(i) ~=0
        cnt=cnt+1;
        loc = mod(PDR.TEST_loc_sequence(i),5);
        scale = PDR.TEST_scale_sequence(i);
        %r= scale;
        sz(i)=ceil(log10(scale)+1)+5;
        colr = find(scale_list==PDR.TEST_scale_sequence(i))/length(scale_list);
        plot(x,PDR.TEST_loc_sequence(i),...
            'Marker','*',...
            'MarkerFaceColor','none',...
            'MarkerEdgeColor',[colr 0 1-colr],...
            'MarkerSize',sz(i));
    end
    
end

PDR.n_test_trials = cnt;

xlabel('Trial #');
ylabel('Location (el,az)');
y_tick_labels{1,PDR.TEST_nlocs} = [];
for i=1:PDR.TEST_nlocs
    y_tick_labels{i} = ['(' num2str(PDR.TEST_locs(i,1)) ',' num2str(PDR.TEST_locs(i,2)) ')'];
end
set(gca,'YTick',1:PDR.TEST_nlocs,'YTickLabel',y_tick_labels);
axis([0 PDR.ntrials+1 0 PDR.TEST_nlocs+1]);
whitebg(gcf,'k');