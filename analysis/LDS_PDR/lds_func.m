function [varargout] = lds_func(fcn,varargin)
global PDR passmein A

switch nargin
    case 4
        var1=varargin{1};
        var2=varargin{2};
        var3=varargin{3};
        feval(fcn,var1,var2,var3);
    case 2
        var1=varargin{1};
        feval(fcn,var1);
    otherwise
        feval(fcn);
end

function load_and_filter_data
    global PDR passmein A
    A.numpts=passmein.buf_pts/(2^passmein.decimateFactor);
    A.numbufs=passmein.nptsTotalPlay/(2^passmein.decimateFactor)/A.numpts;  
    button = questdlg('Use REC1 or REC2???','Choose recording filename','REC1','REC2','REC1');
    fid=fopen([A.pname PDR.filename '_' button A.ext]); A.raw_data=fread(fid,'short'); fclose(fid);
    A.trialtype=A.raw_data(1:A.numpts:size(A.raw_data));
    A.idx2trials=find(A.trialtype==1);
    A.idx2trials = A.idx2trials(1:end); 
    PDR.ntrials = length(A.idx2trials);
    PDRdata.trace=cell(length(A.idx2trials),1);
    % remove last trial if there's a buffer mismatch
    while size(A.raw_data,1) < ((A.numpts*(A.idx2trials(end)-1)+1) + A.bufs_after*A.numpts)
        A.idx2trials=A.idx2trials(1:(end-1));
    end
    len=length(A.idx2trials);
    TRIALS = nan(len,A.bufs_total*(A.numpts-A.infopts)); % temporary matrix for parsing trials

    t=0;
    while t < len
        t=t+1;
        q=1 + A.numpts*(A.idx2trials(t)-1);
        q=q - A.numpts*A.bufs_before; % include A.bufs_before buffers before trial buffer (AC 6/27/2012)
        for k=1:A.bufs_total
            startPt = 1+((k-1)*(A.numpts-A.infopts));
            stopPt = (k*(A.numpts-A.infopts));
            if (q+1) < 1
                TRIALS(t,startPt:stopPt) = nan(A.numpts-A.infopts,1);
            else
                TRIALS(t,startPt:stopPt) = A.raw_data(q+A.infopts:q+A.numpts-1);
            end
            q=q+A.numpts;
        end
        
        % sets values in "prior" buffers equal to first value in 1st
        % buffer (only for first stimulus trace):
        tmp = isnan(PDRdata.trace{t,1});
        tmp = find(tmp == 1);
        if ~isempty(tmp)
            TRIALS(t,tmp) = TRIALS(t,tmp(end)+1).*ones(1,length(tmp));
        end
        PDRdata.numpts(t)=A.bufs_total*(A.numpts-A.infopts);
    end

    PDRdata.trace = mat2cell(TRIALS,ones(1,len));

    for n = 1:length(A.idx2trials)
        P_data (n,:) = PDRdata.trace{n}(1,:)';
    end

    P_data = P_data';
    A.ntrials = length(A.idx2trials);
    y = zeros(A.bufs_total*(A.dec_pts-A.infopts),A.ntrials)*NaN;
    A.meandata = zeros(A.bufs_total*(A.dec_pts-A.infopts),A.ntrials)*NaN;
    %filter AC noise out of trace; normalize all trial DC levels
    h0 = waitbar(0,'Hang on, filtering data...');
    % 20 Hz low pass filter:
    FilterSpec = fdesign.lowpass('N,Fc,Ap,Ast',300,20,1,60,A.dec_Fs);
    FilterObj = design(FilterSpec,'equiripple');
    b=FilterObj.Numerator;
    for trial = 1:A.ntrials
        %Fs = 1000;
        %t = 0:1/Fs:1;
        %b = ones (1,94)/94;
        y (:,trial)= filtfilt(b,1, P_data(:,trial));
        A.meandata (:, trial) = (y(:, trial))-((((y(round(A.bufs_before*(A.dec_pts-A.infopts)+A.sndStart*(A.dec_pts-A.infopts)), trial))))); % changed to zero at sound onset (A.sndStart)
        P_data (:, trial) = ( P_data(:, trial))-(((( P_data(A.bufs_before*(A.dec_pts-A.infopts)+round(A.sndStart*(A.dec_pts-A.infopts)), trial)))));
        waitbar(trial/A.ntrials,h0);
    end
    close(h0)
    if A.invertFlg{1}(1) == 'y'
        A.meandata = -1*A.meandata; %INVERTING ONLY!
    end
    
function remove_blinks
    global PDR passmein A
    %Plot fig to pick blink threshold levels
    hblink = figure('Name','Remove A.trials with Blinks','NumberTitle','off');
    set(gcf, 'Position', [0.1 0.1 0.4*A.scrn(3) 0.4*A.scrn(4)]);
    plot (A.meandata, 'k'); grid; zoom on; set (gcf, 'pointer', 'crosshair');
    
    % Added this to give a marker of sound onset (AC 6/27/2012):
    A.startInt = A.bufs_before*(A.dec_pts-A.infopts)+round((A.dec_pts-A.infopts)*A.sndStart); % start integrating when sound starts in the buffer
    line([A.startInt A.startInt],[min(min(A.meandata)) max(max(A.meandata))]);
    axis tight
    A.lolimit = input ('Lower cutoff =  ');
    A.hilimit = input ('Upper cutoff =  ');
    close(hblink);

    %FILTER OUT ALL A.trials WITH BLINKS

    for trial = 1:A.ntrials
        loblinks = min(find((A.meandata(:,trial)) <= A.lolimit));
        hiblinks = min(find((A.meandata(:,trial)) >= A.hilimit));
        loblinks = isempty (loblinks);
        hiblinks = isempty (hiblinks);
        loblink (trial) = loblinks;
        hiblink (trial) = hiblinks;

    end

    lonay = find (loblink ==1);
    hinay = find (hiblink ==1);
    A.nay = intersect(lonay, hinay);
    A.nbltrials=(length(A.nay));
    
    %cull out blink-free A.trials

    A.mdata = zeros(A.bufs_total*(A.dec_pts-A.infopts),A.nbltrials)*NaN;

    for n = 1:A.nbltrials;
        A.mdata(:,n) = A.meandata(:,A.nay(n));
    end

function pdr_diagnostic_plots
    global PDR passmein A
    %PDR magnitude

    A.startInt = A.bufs_before*(A.dec_pts-A.infopts)+round((A.dec_pts-A.infopts)*A.sndStart); % start integrating when sound starts in the buffer
    
    for m = 1:A.nbltrials;
            A.areas(m) = sum (A.mdata(round(A.startInt):round(A.startInt+A.time2integrate),m)) ;
    end
    
%     hscatt = figure('Name','PDR magnitude per trial','NumberTitle','off');
%     set(gcf, 'Position', [0.05 0.05 0.25*A.scrn(3) 0.25*A.scrn(4) ]);
    subplot(A.h(1));
    plot(A.areas(1),'r*'); hold on; grid on;
    plot (A.areas, 'k.'); 
    
    axis tight
    xlabel ('Trial number'); ylabel ('PDR magnitude mV.s');
    title(['PDR magnitude per trial']);
    legend({'1st Trial','Other Trials'},'Location','SouthOutside')
%     set(A.hFig, 'InvertHardCopy', 'off');
%     if A.pflag
%         print(gcf,'-dpsc2');
%     end
    
    % was there a PDR response?
%     hbatav = figure ('Name','Batch Averages of PDR traces across session',...
%        'NumberTitle','off');
    subplot(A.h(3));
%     set(gcf, 'Position', [0.3*A.scrn(3) 0.05 .25*A.scrn(3) .25*A.scrn(4) ]);
    hold on;
    batches = 5;
    names = cell(1,batches);

    minVal = 0;
    maxVal = 0;
    
%     coloridx = -(0.5*batches):0.5*batches;
%     for i =-(0.5*batches):0.5*batches
%         redval = (-0.5*batches + i)/-batches;
%         grnval = (batches - 2*abs(i))/batches;
%         bluval = (0.5*batches + i)/batches;
%         colorvals(find(coloridx==i),:) = [redval grnval bluval];
%     end
%     
    background_colval=[9.529411764705882e-01 8.705882352941177e-01 7.333333333333333e-01];
    cnt = 1;
    stylez={'-.','--','-'};
    colz=cool;
    for p = 1:batches
        if cnt > 4
            cnt = 1;
        end

        lo(p) = ((floor(A.nbltrials/batches))*(p-1))+1;
        hi(p) = floor(A.nbltrials/batches)*p;
        %plot (mean(A.mdata(:,lo(p):hi(p))'),...
        %   'Color',colorvals(p,:),'LineWidth',2);
        plot(mean(A.mdata(:,lo(p):hi(p))'),'Color',colz(round((p/batches)*64),:),'LineWidth',2,'LineStyle',stylez{mod(p,3)+1});
        
        minVal = min(minVal,min(mean(A.mdata(:,lo(p):hi(p))')));
        maxVal = max(maxVal,max(mean(A.mdata(:,lo(p):hi(p))')));
        
        names{p} = [num2str(lo(p)) '-' num2str(hi(p))];
        cnt = cnt + 1;
    end

    % Added this to give a marker of sound onset (AC 6/27/2012)
    line([A.startInt A.startInt],[minVal maxVal],'Color','r');
    axis tight

    endtime = round((A.bufs_total*A.bufftime) - A.sndStart*A.bufftime - A.bufs_before*A.bufftime);
    begintime = -round(A.bufs_before*A.bufftime + A.sndStart*A.bufftime);
    set (gca, 'xtick', [0:(A.dec_pts-A.infopts):A.bufs_total*(A.dec_pts-A.infopts)], 'xticklabel',[begintime:round(endtime/A.bufs_total):endtime]);
    set(gca,'Color',background_colval);
    xlabel ('Time after sound onset (ms)'); ylabel ('Pupil size (mV)');
    title (['Batch Average Traces']);
    legend(names,'Location','SouthOutside','Color',background_colval)
    grid off;
    
function running_avg_habtrials(winnum,windur,subp)
    global PDR passmein A
    % winnum is the window (# of A.trials to average per window)
    % windur is the window length in minutes
    A.startInt = A.bufs_before*(A.dec_pts-A.infopts)+round((A.dec_pts-A.infopts)*A.sndStart); % start integrating when sound starts in the buffer
    
    % calc initial windowed avg
    in = 0;
    area_list = zeros(round(winnum),1)*NaN;
    for i = 1:round(winnum)
        in = in + sum(A.mdata(round(A.startInt):round(A.startInt+A.time2integrate),i));
        area_list(i) = sum(A.mdata(round(A.startInt):round(A.startInt+A.time2integrate),i));
    end
    
    % stdev and mean for first window
    begin_std = std(area_list);
    begin_avg = in/round(winnum);

    % running avg for all other windows
    in = 0;
    avg = zeros(A.nbltrials,1)*NaN;
    stdevs = zeros(A.nbltrials,1)*NaN;
    stdevs(1) = begin_std;
    avg(1) = begin_avg;
    for j = 2:round(A.nbltrials-winnum)
        area_list = zeros(round(winnum),1)*NaN;
        for k = 1:round(winnum)
            in = in + sum(A.mdata(round(A.startInt):round(A.startInt+A.time2integrate),j+k));
            area_list(k) = sum(A.mdata(round(A.startInt):round(A.startInt+A.time2integrate),k));
        end
        stdevs(j) = std(area_list);
        avg(j) = in/round(winnum);
        in = 0;
    end
%
    xes = zeros(A.nbltrials,1)*NaN;
    for m = 1:A.nbltrials
        minute_cnt = floor(m*windur/winnum);
        xes(m) = minute_cnt + mod(m,winnum/windur)*windur/winnum;
    end

    %traceplot = figure('Name',['Running Avg.: ' num2str(windur) ' Min. Window'],'NumberTitle','off');
    %set(gcf, 'Position', [0.55*A.scrn(3) 0.05 0.25*A.scrn(3) 0.25*A.scrn(4) ]);
    subplot(subp);
    plot(xes,avg,'k-'); hold on;
%     plot(xes,avg+stdevs,'k--','LineWidth',2);
%     plot(xes,avg-stdevs,'k--','LineWidth',2);
    xlabel ('Session Time (Min)'); ylabel ('Avg Response');
    title (['Running Avg. Hab Response: Win. Length = ' num2str(windur) ' Min']);
    plotname = ['RunningAvg_' num2str(windur) 'MinWindow'];

function parse_trials
    global PDR passmein A
   
    
    
    % sort A.trials into hab, and each test sound type
    snd = PDR.trialID(1:end)';
    
    
    for n = 1:A.nbltrials;
        nblsnd(n) = snd(A.nay(n));
    end

    %get list of trial types (hab and other)
    A.trials.types = unique(nblsnd);
    A.trials.types = A.trials.types(find(A.trials.types~=A.habID));
    A.trials.types = sort(A.trials.types);

    %get data for hab A.trials
    A.trials.hab.idx = find(nblsnd==A.habID);
    A.trials.hab.areas = A.areas(A.trials.hab.idx);
    A.trials.hab.traces = A.mdata(:,A.trials.hab.idx);

    %get overall mean and std deviation
    A.trials.stdev = std(A.areas);
    A.trials.mean = mean(A.areas);

    %get zscores for hab A.trials
    A.trials.hab.zscores = (A.trials.hab.areas - A.trials.mean)/A.trials.stdev;

    for i = 1:length(A.trials.types)
        A.trials.test{i}.loc = [num2str(PDR.SOUNDS_azimuths(A.trials.types(i)))];
        A.trials.test{i}.idx = find(nblsnd==A.trials.types(i));
        A.trials.test{i}.areas = A.areas(A.trials.test{i}.idx);
        A.trials.test{i}.traces = A.mdata(:,A.trials.test{i}.idx);
        A.trials.test{i}.zscores = (A.trials.test{i}.areas - A.trials.mean)/A.trials.stdev;
    end
    
function plot_average_traces
    global PDR passmein A
 
    A.startInt = A.bufs_before*(A.dec_pts-A.infopts)+round((A.dec_pts-A.infopts)*A.sndStart); % start integrating when sound starts in the buffer
    
    %traceplot = figure('Name','Averaged Test and Hab A.trials','NumberTitle','off');
    %set(gcf, 'Position', [0.05 0.3*A.scrn(4) .25*A.scrn(3) .25*A.scrn(4) ]);
    grid off
    plot(mean(A.trials.hab.traces'),'w','LineWidth',2);
    minVal = min(min(mean(A.trials.hab.traces')));
    maxVal = max(max(mean(A.trials.hab.traces')));
    
    hold on
    leg = {'Hab Trials'};
    back_colr = [.5 .5 .5]; %[7.294117808341980e-01 8.313725590705872e-01 9.568627476692200e-01];
    set(gca,'Color',back_colr);
%     plot(mean(A.trials.hab.traces(:,100:end)'),'k-','LineWidth',2);
%     leg{2} = 'Hab A.trials (w/o 1st 100)';
%     
%     plot(mean(A.trials.hab.traces(:,ceil(end/2):end)'),'k--','LineWidth',2);
%     leg{3} = 'Hab A.trials (Last Half)';

    tt_types = (length(A.trials.types));
    stylez={'-','--','-.'};
    colz=hsv;
    for i = 1:tt_types
        plot(mean(A.trials.test{i}.traces'),'Color',colz(round((i/tt_types)*64),:),'LineWidth',2,'LineStyle',stylez{mod(i,3)+1});
        leg{i+1} = ['\Delta = ' num2str(abs(PDR.LAG_hab_pos - str2num(A.trials.test{i}.loc))) '^{\circ}'];
        minVal = min(minVal,min(min(mean(A.trials.test{i}.traces'))));
        maxVal = max(maxVal,max(max(mean(A.trials.test{i}.traces'))));
    end
    
    legend(leg,'Location','EastOutside','Color',back_colr);
    
    % Added this to give a marker of sound onset (AC 6/27/2012):
    len = size(A.trials.hab.traces,1);
    line([A.startInt A.startInt],[minVal maxVal],'Color','r');
    axis tight
    endtime = (A.bufs_total*A.bufftime) - A.sndStart*A.bufftime - A.bufs_before*A.bufftime;
    begintime = -(A.bufs_before*A.bufftime + A.sndStart*A.bufftime);
    
    for k=0:A.bufs_total
        labl{k+1}=num2str(begintime+A.bufftime*k,'%.1f');
    end
    set(gca, 'XTick', [0:(A.dec_pts-A.infopts):A.bufs_total*(A.dec_pts-A.infopts)], 'XTickLabel',labl);

    xlabel ('Time after sound onset (ms)'); ylabel ('Pupil size (mV)');
    title (['Average traces of Habituating and Test Trials']);
    
function plot_trial_traces
    global PDR passmein A
    
    A.startInt = A.bufs_before*(A.dec_pts-A.infopts)+round((A.dec_pts-A.infopts)*A.sndStart); % start integrating when sound starts in the buffer


    
    tt_types=(length(A.trials.types)-1);
    for i = 1:tt_types
        
        leg = cell(1,1);
        traceplot = figure('Name',['Individual A.trials: ' A.trials.test{i}.loc ' Correlation'],'NumberTitle','off');
        set(gcf, 'Position', [.3*A.scrn(3) 0.3*A.scrn(4) .25*A.scrn(3) .25*A.scrn(4) ]);
        subplot(A.h(1));
        set(gca,'Color',[.8 .8 .8]);
        plot(mean(A.trials.hab.traces'),'k:','LineWidth',2);
        
        minVal = min(min(mean(A.trials.hab.traces')));
        maxVal = max(max(mean(A.trials.hab.traces')));
        
        hold on
        leg{1} = 'Hab A.trials';
        plot(mean(A.trials.hab.traces(:,100:end)'),'k-','LineWidth',2);
        leg{2} = 'Hab A.trials (w/o 1st 100)';
        plot(mean(A.trials.hab.traces(:,ceil(end/2):end)'),'k--','LineWidth',2);
        leg{3} = 'Hab A.trials (Last Half)';
        
        
        num_trials = size(A.trials.test{i}.traces,2);
        prior_hab_trials_1MinWin_traces = [];
        prior_hab_trials_5MinWin_traces = [];

        
        for j = 1:num_trials
            
            ind = find(A.trials.hab.idx < A.trials.test{i}.idx(j));
            prior_hab_trials_1MinWin_idx = find(A.trials.hab.idx(ind) ...
                > (A.trials.test{i}.idx(j) - A.win1) ...
                ); % find hab A.trials during 1 minute window prior to this test trial
            prior_hab_trials_5MinWin_idx = find(A.trials.hab.idx(ind) ...
                > (A.trials.test{i}.idx(j) - A.win5) ...
                ); % find hab A.trials during 5 minute window prior to this test trial

            prior_hab_trials_1MinWin_traces = ...
                [prior_hab_trials_1MinWin_traces A.mdata(:,prior_hab_trials_1MinWin_idx)]; % get traces for hab A.trials 1 minute prior to this test trial
            prior_hab_trials_5MinWin_traces = ...
                [prior_hab_trials_5MinWin_traces A.mdata(:,prior_hab_trials_5MinWin_idx)]; % get traces for hab A.trials 5 minutes prior to this test trial
            
            plot(A.trials.test{i}.traces(:,j),'Color',[1 j/num_trials j/num_trials],'LineWidth',3);
            minVal = min(minVal,min(min(A.trials.test{i}.traces(:,j))));
            maxVal = max(maxVal,max(max(A.trials.test{i}.traces(:,j))));
            leg{j+3} = ['Trial #' num2str(j)];
            cnt = j+3;
        end
        
        leg{cnt+1} = 'Hab A.trials (1 Min Prior to Test Trial)';
        plot(mean(prior_hab_trials_1MinWin_traces'),'m-','LineWidth',2);
        leg{cnt+2} = 'Hab A.trials (5 Min Prior to Test Trial)';
        plot(mean(prior_hab_trials_5MinWin_traces'),'m--','LineWidth',2);
        
        %set(gca,'Color',[.7 .7 .7]);
        legend(leg);
        grid on;
        
        % Added this to give a marker of sound onset (AC 6/27/2012):
        len = size(A.trials.hab.traces,1);
        line([A.startInt A.startInt],[minVal maxVal],'Color','r');
        axis tight
        endtime = round((A.bufs_total*A.bufftime) - A.sndStart*A.bufftime - A.bufs_before*A.bufftime);
        begintime = -round(A.bufs_before*A.bufftime + A.sndStart*A.bufftime);
        set (gca, 'xtick', [0:(A.dec_pts-A.infopts):A.bufs_total*(A.dec_pts-A.infopts)], 'xticklabel',[begintime:round(endtime/A.bufs_total):endtime]);
        
        
        xlabel ('Time after sound onset (ms)'); ylabel ('Pupil size (mV)');
        title (['Individual ' A.trials.test{i}.loc ' Correlation Test']);
        
        plotname = ['Traces_' A.trials.test{i}.corr(1:end-1) '_PercentCorr'];
    end

function plot_zscores
        
        global A PDR passmein
        
        %traceplot = figure('Name','Plot of z-scores for Each Condition','NumberTitle','off');
        %set(gcf, 'Position', [0.15 0.15 0.4*A.scrn(3) 0.4*A.scrn(4) ]);
        grid on;
        len = 0;
        for i = 1:(length(A.trials.types)-1)
            len = len + size(A.trials.test{i}.zscores,2);
        end
        
        len = len + size(A.trials.hab.zscores,2);
        xes = zeros(len,1)*NaN;
        zes = zeros(len,1)*NaN;
        
        box_labels = cell(len,1);
        xes(1:size(A.trials.hab.zscores,2)) = PDR.LAG_hab_pos*ones(size(A.trials.hab.zscores,2),1);
        len2 = size(A.trials.hab.zscores,2)+1;
        for p=1:size(A.trials.hab.zscores,2)
            pos = num2str(PDR.LAG_hab_pos);
            box_labels{p} = pos;
        end
        zes(1:size(A.trials.hab.zscores,2)) = A.trials.hab.zscores;
        
        for j = 1:length(A.trials.types)
            xes(len2:(len2+size(A.trials.test{j}.zscores,2)-1)) = str2num(A.trials.test{j}.loc)*ones(size(A.trials.test{j}.zscores,2),1);
            zes(len2:(len2+size(A.trials.test{j}.zscores,2)-1)) = A.trials.test{j}.zscores;
            pos = num2str(A.trials.test{j}.loc);
            for p = len2:(len2+size(A.trials.test{j}.zscores,2)-1)
                box_labels{p} = pos;
            end
            len2 = len2 + size(A.trials.test{j}.zscores,2);
            
        end
%         boxplot(zes,box_labels);
%         figure;
        hab_end=size(A.trials.hab.zscores,2);
        %x_jitter=-1+2.*rand(length(xes(1:hab_end)),1);
        plot(xes(1:hab_end),zes(1:hab_end),'k*'); hold on;
        %x_jitter=-1+2.*rand(length(xes(hab_end+1:end)),1);
        plot(xes(hab_end+1:end),zes(hab_end+1:end),'m+','MarkerSize',7); %'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 0 0]);
        axis([min(xes)-5 max(xes)+5 min(zes(hab_end+1:end))-1 max(zes(hab_end+1:end))+1]);
%         hold on;
%         
%         begin_pt = size(A.trials.hab.zscores,2)+1; % not using hab trials for regression fit
%         p = polyfit(xes(begin_pt:end),zes(begin_pt:end),1);
%         zfit = polyval(p,xes(begin_pt:end));
%         plot(xes(begin_pt:end),zfit,'r-');
%         
%         resids = zes(begin_pt:end) - zfit;
%         SSresid = sum(resids.^2);
%         SStotal = (length(zes(begin_pt:end))-1)*var(zes(begin_pt:end));
%         rsq = 1 - SSresid/SStotal;
%         
%         set(gca,'Color',[.8 .8 .8]);
%         grid on; axis tight;
        set(gca,'Color',[.7 .7 .7]);
        xlabel ('Azimuth'); ylabel ('z-scores');
         title (['z-scores for each location tested']);

function roc_PERFCURVE(i)

global PDR A passmein

%NOTE: MUST USE MATLAB VERSION R2010a OR LATER FOR PERFCURVE() TO WORK

% ANDREW (Note): Takes a list of pdr response areas
% ("areahabit" for habituated probes and "areatest" for test probes) for a
% particular test condition and computes percent correct.
% then plots roc function for the condition

% function [percent_correct] = roc_percent(areahabit,areatest,rocfig,corr_val,total_tests);
% gets percent area under the ROC curve, given areahabit and areatest

% rocfig = handle to roc plot figure
% corr_val = correlation value of a particular test condition
% total_tests = total # of test conditions
warning off stats:perfcurve:SubSampleWithMissingClasses
testID=A.test_info{i};
colr=[1 0 0]; %[i/length(A.trials.types) 1 (1-(i/length(A.trials.types)))];
%mkrs={'s','d','o'};
%stylez={'-','-.','--'};
%style_type=stylez{mod(i,3)+1};
%mkr=mkrs{mod(i,3)+1};
areatest=A.trials.test{i}.areas;
areahabit = randsample(A.trials.hab.areas,100);
scores = [areahabit areatest]; % combined areas

labels(1:length(areahabit)) = 0;
labels(length(areahabit)+1:length(scores)) = 1;

[X,Y,T,AUC] = perfcurve(labels,scores,1,'NBoot',10);

figure(A.rocfig); hold on;
hSub=subplot(ceil(A.total_tests/A.bufs_before),4,i);
%X2=[X(:,1)'; X(:,1)'; X(:,1)'];
%Y2=[Y(:,2)' - Y(:,1)'; Y(:,1)'; Y(:,3)' - Y(:,1)'];
h=area(Y(:,1)'); set(h,'FaceColor',[i/length(A.trials.types) 0 0]);
%colormap summer
%set(h,'FaceColor',[i/length(A.trials.types) 0 0]);
%set(h,'Marker','none','Color',colr,'LineWidth',3); %'LineStyle',style_type);
% plot errors
%errorbar(X(:,1),Y(:,1),Y(:,2)-Y(:,1),Y(:,3)-Y(:,1));
axis([1 length(Y(:,1)') 0 1]);
len=length(Y(:,1)');
set(hSub,'XTick',1:(len/2):len,'XTickLabel',{'0','0.5','1'});
xlabel('False positive rate'); ylabel('True positive rate');
A.pc(i) = AUC(1);
A.lowerCI(i) = AUC(2);
A.upperCI(i) = AUC(3);
title(['Test Separation = ' num2str(abs(PDR.LAG_hab_pos - str2num(A.trials.test{i}.loc))) ' Azimuth PC=' num2str(AUC(1))]);
drawnow;

function running_roc
global PDR A passmein

%NOTE: MUST USE MATLAB VERSION R2010a OR LATER FOR PERFCURVE() TO WORK

N=round(A.dec_Fs);
D=10; % resamples to ~10Hz samp. rate
step=round(N/D);
% determine normalization factor at each time point
nTestTypes=length(A.test_info);
tmp = A.trials.hab.traces;
for j0=1:nTestTypes
    tmp = [tmp A.trials.test{j0}.traces];
end
norm_factor = std(tmp,0,2); % stdev at each time point

ROC = nan(nTestTypes,length((A.startInt+1):step:(A.startInt+A.time2integrate)));
LowerCI=ROC; UpperCI=ROC;
hWait=waitbar(0,'crunching...'); cnt=0;

for tm=(A.startInt+1):step:(A.startInt+A.time2integrate)
    clear scores_habit
    cnt=cnt+1;
    scores_habit = randsample(A.trials.hab.traces(tm,:),100);
    for k=1:nTestTypes
        clear scores_test scores labels
        scores_test=A.trials.test{k}.traces(tm,:);
        scores = [scores_habit scores_test];
        labels = [zeros(1,length(scores_habit)) ones(1,length(scores_test))];
        [X,Y,T,AUC] = perfcurve(labels,scores,1);
        ROC(k,cnt)=AUC(1);
        %LowerCI(k,cnt)=AUC(2);
        %UpperCI(k,cnt)=AUC(3);
    end
    waitbar(cnt/length((A.startInt+1):step:(A.startInt+A.time2integrate)));
end
close(hWait);
figure;
errorbar(mean(ROC,2),std(ROC,0,2));
% dim = ceil(sqrt(nTestTypes));
% alpha=0.5; colz=hsv;
% F=(1:length((A.startInt+1):step:(A.startInt+A.time2integrate)))./A.dec_Fs;
% for j2=1:nTestTypes
%     subplot(dim,dim,j2);
%     acolor = colz(round((j2/nTestTypes)*64),:);
%     fill([F fliplr(F)],[LowerCI(j2,:) fliplr(UpperCI(j2,:))],acolor,'linestyle','none','FaceAlpha', alpha);
%     hold on;
%     plot(F,ROC(j2,:),'Color',acolor,'linewidth',1.5);
% end

function plot_stdshade_avg_traces
    global PDR passmein A
    
    begintime = -(A.bufs_before*A.bufftime + A.sndStart*A.bufftime);
    endtime = (A.bufs_total*A.bufftime) - abs(begintime);
    
    for k=0:A.bufs_total
        labl{k+1}=num2str(begintime+A.bufftime*k,'%.1f');
    end
    
    back_colr = [.5 .5 .5]; %[7.294117808341980e-01 8.313725590705872e-01 9.568627476692200e-01];
    alpha = 0.5;
    tt_types = (length(A.trials.types));
    colz=hsv;
    dim = ceil(sqrt(tt_types));
    minVal = min(min(mean(A.trials.hab.traces') - std(A.trials.hab.traces')));
    maxVal = max(max(mean(A.trials.hab.traces') + std(A.trials.hab.traces')));
    for i0 = 1:tt_types
        subplot(dim,dim,i0); hold on;
        title(['\Delta = ' num2str(abs(PDR.LAG_hab_pos - str2num(A.trials.test{i0}.loc))) '^{\circ}']);
        if length(A.trials.test{i0}.zscores)>1
            amatrix = A.trials.hab.traces';
            acolor='w';
            stdshade(amatrix,alpha,acolor);
            hold on
            set(gca,'Color',back_colr);
            amatrix = A.trials.test{i0}.traces';
            acolor = colz(round((i0/tt_types)*64),:);
            stdshade(amatrix,alpha,acolor);
            minVal = min(minVal,min(min(mean(A.trials.test{i0}.traces') - std(A.trials.test{i0}.traces'))));
            maxVal = max(maxVal,max(max(mean(A.trials.test{i0}.traces') + std(A.trials.test{i0}.traces'))));
            set(gca, 'XTick', [0:(A.dec_pts-A.infopts):A.bufs_total*(A.dec_pts-A.infopts)], 'XTickLabel',labl);
        end
    end
    
    for i1=1:tt_types
        if length(A.trials.test{i0}.zscores)>1
        subplot(dim,dim,i1); hold on;
        line([A.startInt A.startInt],[minVal maxVal],'Color','w');
        axis([0 A.bufs_total*(A.dec_pts-A.infopts) minVal maxVal]);
        xlabel ('Time after sound onset (s)'); ylabel ('Pupil size (mV)');
        end
    end
