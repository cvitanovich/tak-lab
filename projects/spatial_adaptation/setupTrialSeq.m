function quit = setupTrialSeq()
% function to set up trial sequence for spatialPDR
global PDR session

setup=1;
while setup
	% Setup randomized scale and loc sequences:
	scales = PDR.TEST_scales;
	
	mu = mean(scales);
	scales = sort(scales);
	
	% lowest test scales:
	lows = scales(1:round(length(scales)*.333)); % lower third
	scale_seq = zeros(1,PDR.ntrials+1);
	start = PDR.npretrials+1;
	stop = PDR.ntrials+1 - mod(PDR.ntrials-PDR.npretrials+1,PDR.TEST_trial_freq);
	nm = floor((stop - start)/PDR.TEST_trial_freq)+1;
	tmp=ceil(length(lows)*rand(1,nm));
	scale_seq(start:PDR.TEST_trial_freq:stop) = lows(tmp);
	
	% middle scales:
	mids = scales(round(length(scales)*.333)+1:round(length(scales)*.667)); % middle third
	% highest scales:
	highs = scales(round(length(scales)*.667)+1:end);
	% combine mids and highs
	mid_hi = [mids highs];
	
	% randomized sequence with this pattern:
	% (hi or med) low (hi or med) low (hi or med) ...
	FRQ = 2*PDR.TEST_trial_freq;
	start = PDR.npretrials+1;
	stop = PDR.ntrials+1 - mod(PDR.ntrials-PDR.npretrials+1,FRQ);
	nm = floor((stop - start)/FRQ)+1;
	tmp = ceil(length(mid_hi)*rand(1,nm));
	scale_seq(start:FRQ:stop) = mid_hi(tmp);
	
	PDR.TEST_scale_sequence = scale_seq;
	
	% create a sequence of locations that is randomized:
	start = PDR.npretrials+1;
	stop = PDR.ntrials+1 - mod(PDR.ntrials-PDR.npretrials+1,PDR.TEST_trial_freq);
	nm = floor((stop - start)/PDR.TEST_trial_freq)+1;
	
	% randomized sequence of locations:
	loc_seq = zeros(1,PDR.ntrials+1);
	for i = start:PDR.TEST_trial_freq:stop
        loc_seq(i) = round(1 + (PDR.TEST_nlocs - 1) * rand(1)); % gives a random integer indicating trial loc IDs
	end
	
	PDR.TEST_loc_sequence = loc_seq(1:PDR.ntrials+1);


    screen_size = get(0, 'ScreenSize');
	session.hTrialPlot = figure;
	set(session.hTrialPlot, 'Position', [0.1*screen_size(3) 0.1*screen_size(4) 0.7*screen_size(3) 0.7*screen_size(4)] );
	hold on;
	x=0;
	tmp = find(PDR.TEST_loc_sequence~=0);
    if length(tmp)<1
        warndlg('No trials with this config!','No Trials');
        quit = 1;
        close(session.hTrialPlot);
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
            sz(i)=ceil(log10(scale)+1)*2+4;
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