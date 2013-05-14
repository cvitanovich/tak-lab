% lds daily analysis
clear all;
global A
A.pflag = 0; %print flag (set to 1 to print figures)
A.sflag = 1; %save flag (set to 1 to save figures)
A.ext = '.frf';

% SELECT HEADER FILE FOR ANALYSIS

current_dir = pwd;
dt='/home/andrew/science/Alex/LDS_Project/analysis/code';
cd ('/home/andrew/science/Alex/LDS_Project/data');
path(path,dt);
[A.fname, A.pname] = uigetfile('*.mat', 'Select header file for daily analysis');
global PDR passmein

if A.fname==0
    disp('Couldn''t find file!')
    return
end

A.analysis_dir = ['/home/andrew/science/Alex/LDS_Project/analysis/data/'];

cd (A.pname)
load(A.fname)
if ~isempty(dt)
    A.DATApath = dt;
end
for tID=1:length(PDR.TEST_azimuths)
    A.test_info{tID} = num2str(PDR.TEST_azimuths(tID));
end
A.infopts=2; % first 2 points used for trial info (trial marker=1, location=1,2,3...,8,9,etc)
A.PDRfile=passmein.FN1;
tmp = find(PDR.LAG_sounds{1}~=0);
A.sndStart = tmp(1)/size(PDR.LAG_sounds{1},2); % how far through the buffer (fractional) does the sound start
A.dec_Fs = PDR.stim_Fs/2^PDR.decimationfactor; % decimated samp. rate
A.dec_pts = PDR.buf_pts/2^PDR.decimationfactor; % num buf pts (decimated)
A.time2integrate = A.dec_Fs*2; % using integration time of 2 sec
A.bufftime = A.dec_pts/A.dec_Fs; % actual timelength of buffer in seconds
A.defaultFlg = {'n'};
A.invertFlg=inputdlg('Invert traces? (y/n)','Invert Trace Dlg',1,A.defaultFlg);
A.bufs_before=4; A.bufs_after=4;
A.bufs_total=A.bufs_before+A.bufs_after;
A.scrn = get(0, 'ScreenSize');
A.file='';
for i=1:length(A.PDRfile)
    if A.PDRfile(i)~='_'
        A.file=[A.file A.PDRfile(i)];
    else
        A.file=[A.file '\' A.PDRfile(i)];
    end
end

%***************************************************
%*******LOAD AND FILTER DATA ***********************
%***************************************************

lds_func('load_and_filter_data');

if length(A.idx2trials) > PDR.ntrials
    PDR.ntrials = length(A.idx2trials);
end


%***************************************************
%*******REMOVE BLINKS ******************************
%***************************************************
lds_func('remove_blinks');

if A.sflag
    cd (A.analysis_dir) % move to analysis dir to save figs, data, etc.
end


% MAIN FIGURE SETUP

A.hFig=figure;
set(gcf,'Color',[1 1 1],'Position',[.02*A.scrn(3) .04*A.scrn(4) .20*A.scrn(3) .9*A.scrn(4)]);
rowz=11; colz=8;
% trial mag subplot
A.h(1) = subplot(rowz,colz,[73:76 81:84]);
% avg trace each test condition
A.h(2) = subplot(rowz,colz,[1:24]);
% batch avg plots
A.h(3) = subplot(rowz,colz,[37:40 45:48 53:56 61:64]);
% zscores
A.h(4) = subplot(rowz,colz,[33:35 41:43 49:51 57:59]);
% roc summary
A.h(5) = subplot(rowz,colz,[78:80 86:88]);

%***************************************************
%*******DIAGNOSTIC PLOTS ***************************
%***************************************************

% Plot PDR magnitudes for all trials
% Plot Batch averages of hab trials

lds_func('pdr_diagnostic_plots');

%*****************************************************
% parse trials (and sort zscores, areas, traces, etc.)
%*****************************************************
A.habID=find(PDR.SOUNDS_azimuths==PDR.LAG_hab_pos);
for i0=1:length(PDR.SOUNDS_location_sequence)
    PDR.trialID(i0)=find(PDR.SOUNDS_location_sequence(2,i0)==PDR.SOUNDS_azimuths);
end
if length(PDR.SOUNDS_location_sequence)<PDR.ntrials
    for i1=(length(PDR.SOUNDS_location_sequence)+1):PDR.ntrials
        PDR.trialID(i1)=A.habID; % this is a hab trial
    end
end
lds_func('parse_trials');



%***************************************************
% plot average traces
%***************************************************
subplot(A.h(2));
lds_func('plot_average_traces');
A.hShade=figure;
set(gcf,'Color',[1 1 1],'Position',[.6*A.scrn(3) .04*A.scrn(4) .40*A.scrn(3) .9*A.scrn(4)]);
lds_func('plot_stdshade_avg_traces');
PDR
figure(A.hFig);

if length(A.trials.types) == 1 % there were only habit trials
    
% do no more analysis
else % there were test trials
    
    %***************************************************
    % plot zscores and fit a line
    %***************************************************
    subplot(A.h(4));
    lds_func('plot_zscores');
    grid on;


    %***************************************************
    % roc analysis
    %***************************************************
    A.rocfig= figure('Name','ROC Plots','NumberTitle','off');
    set(A.rocfig, 'Position', [0.15 0.15 0.4*A.scrn(3) 0.8*A.scrn(4) ]);
    hold on;
    A.loc=NaN*ones(length(A.trials.types),1);
    A.total_tests=length(A.trials.types);
    for i = 1:length(A.trials.types)
        lds_func('roc_PERFCURVE',i);
        A.loc(i) = str2num(A.trials.test{i}.loc);
    end
    %whitebg(gcf,[.5 .5 .5])
    %legend(A.test_info,'Location','Best')
    
    set(A.rocfig, 'InvertHardCopy', 'off');
%     if pflag
%         print(gcf,'-dpsc2');
%     end

    %hROC_summary = figure('Name','ROC Summary','NumberTitle','off'); 
    %set(gcf, 'Position', [0 0 A.scrn(3) A.scrn(4) ]);
    figure(A.hFig); subplot(A.h(5));
    plot(A.loc,A.pc,'s','MarkerSize',12,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
    hold on;
    %errorbar(corrs,pc,upperCI-pc,lowerCI-pc); DON'T PLOT ERROR BARS
    % TOO FEW TRIALS TO PLOT ERROR BARS

    xlabel('Locations'); ylabel('Percent Correct');
    %axis([-5 100 0 1]);
    title(['Summary Plot of ROC Analysis']);
    
    
end


title_txt{1}=['Analysis for session: ' A.file];
annotation(A.hFig,'textbox',...
    [0.10 0.95 0.8 0.03],...
    'String',title_txt,...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Lucida',...
    'FitBoxToText','off');



% save analysis summary figure and analysis data in analysis dir
if A.sflag
    cd(A.analysis_dir);
    saveas(A.hFig,[A.fname(1:end-4) '_analysis.fig']);
    saveas(A.hFig,[A.fname(1:end-4) '_analysis.png']);
end
    
cd(current_dir);