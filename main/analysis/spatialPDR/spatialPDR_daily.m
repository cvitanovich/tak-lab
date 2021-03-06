% spatialPDR daily analysis
clear all;
global A
A.pflag = 0; %print flag (set to 1 to print figures)
A.sflag = 1; %save flag (set to 1 to save figures)
A.infopts = 1;

button = questdlg('Virtual or Freefield? ','Select Expt Type','Virtual','Freefield','Virtual');
if(strcmp(button,'Virtual'))
    A.ext = '.vrt';
else
    A.ext = '.frf';
end

% SELECT HEADER FILE FOR ANALYSIS

current_dir = pwd;
compID=[];
try
    compID=getenv('USER');
end

% get hostname of system
[~, hname] = system('hostname');

if(strcmp('cvit-macbook',hname(1:12)))
    code_dir='/home/andrew/code/tak-lab/main/analysis/spatialPDR/';
    data_dir='/home/andrew/data/';
    analysis_dir='/home/andrew/analysis_data/';
% elseif(strcmp('cvitanovich',compID)) % my macbook
%     code_dir='/Users/cvitanovich/Documents/MATLAB/tak-lab/analysis/LDS_PDR/';
%     data_dir='/Users/cvitanovich/Desktop/alex_data/';
%     analysis_dir='/Users/cvitanovich/Documents/MATLAB/analysis_data/';
% elseif(strcmp('andrew',compID)) % linux desktop (Arch)
%     %code_dir='/home/andrew/science/Alex/LDS_Project/analysis/code';
%     data_dir='/home/andrew/science/data/projects/AdaptProject/';
%     analysis_dir='/home/andrew/science/analysis/projects/AdaptProject/';
else
    disp('figure out your paths!')
end

[A.fname, A.pname] = uigetfile([data_dir '*.mat'], 'Select header file for daily analysis');
global PDR

if A.fname==0
    disp('Couldn''t find file!')
    return
end

A.analysis_dir = analysis_dir;

%cd (A.pname)
load([A.pname A.fname]);


A.SCALE_list=[PDR.TEST_scales PDR.TEST_outlier_scales];
A.SPL_list=[PDR.TEST_SPLs PDR.TEST_outlier_SPLs];
tmp=find(PDR.TEST_scale_sequence~=0);
A.TEST_level_sequence=-Inf.*ones(1,length(tmp));
for(j=1:length(tmp))
   % find matching SPL for each scale
   A.TEST_level_sequence(tmp(j)) = A.SPL_list( find( A.SCALE_list==PDR.TEST_scale_sequence(tmp(j)) ) );
end

% if ~isempty(code_dir)
%     A.DATApath = code_dir;
% end
for tID=1:length(A.SPL_list)
    A.test_info{tID} = num2str(A.SPL_list(tID));
end
A.PDRfile=PDR.filename;
A.sndStart = PDR.TEST_start_pt/PDR.buf_pts;
%TEST_on_delay_pts/size(PDR.LAG_sounds{1},2); % how far through the buffer (fractional) does the sound start
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

spatialPDR_func('load_and_filter_data');

if length(A.idx2trials) > PDR.ntrials
    PDR.ntrials = length(A.idx2trials);
end


%***************************************************
%*******REMOVE BLINKS ******************************
%***************************************************
spatialPDR_func('remove_blinks');

if A.sflag
    cd (A.analysis_dir) % move to analysis dir to save figs, data, etc.
end

% MAIN FIGURE SETUP

A.hFig=figure;
if(strcmp(compID,'andrew'))
    set(gcf,'Color',[1 1 1],'Position',[.02*A.scrn(3) .04*A.scrn(4) .20*A.scrn(3) .9*A.scrn(4)]);
elseif(strcmp(compID,'cvitanovich'))
    set(gcf,'Color',[1 1 1],'Position',[.02*A.scrn(3) .04*A.scrn(4) .5*A.scrn(3) .9*A.scrn(4)]);
else
    disp('figure out your paths and figure sizes!')
end
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

spatialPDR_func('pdr_diagnostic_plots');

%*****************************************************
% parse trials (and sort zscores, areas, traces, etc.)
%*****************************************************
A.habID=-Inf;

for i0=1:length(PDR.TEST_scale_sequence)
    tmp=find(PDR.TEST_scale_sequence(i0)==A.SCALE_list);
    if(isempty(tmp))
        PDR.trialID(i0)=A.habID;
    else
        PDR.trialID(i0)=A.SPL_list(tmp(1));
    end
end
% if length(PDR.TEST_scale_sequence)<PDR.ntrials
%     for i1=(length(PDR.SOUNDS_location_sequence)+1):PDR.ntrials
%         PDR.trialID(i1)=A.habID; % this is a hab trial
%     end
% end
spatialPDR_func('parse_trials');



%***************************************************
% plot average traces
%***************************************************
subplot(A.h(2));
spatialPDR_func('plot_average_traces');
%A.hShade=figure;
%set(gcf,'Color',[1 1 1],'Position',[.6*A.scrn(3) .04*A.scrn(4) .40*A.scrn(3) .9*A.scrn(4)]);
%spatialPDR_func('plot_stdshade_avg_traces');

figure(A.hFig);

if length(A.trials.types) == 1 % there were only habit trials
    
% do no more analysis
else % there were test trials
    
    %***************************************************
    % plot zscores and fit a line
    %***************************************************
    subplot(A.h(4));
    spatialPDR_func('plot_zscores');
    grid on;

    if(1)
    %***************************************************
    % roc analysis
    %***************************************************
    A.rocfig= figure('Name','ROC Plots','NumberTitle','off');
    set(A.rocfig, 'Position', [0.15 0.15 0.4*A.scrn(3) 0.8*A.scrn(4) ]);
    hold on;
    A.spls=NaN*ones(length(A.trials.types),1);
    A.total_tests=length(A.trials.types);
    for i = 1:length(A.trials.types)
        spatialPDR_func('roc_PERFCURVE',i);
        A.spls(i) = str2num(A.trials.test{i}.spl);
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
    plot(A.spls,A.pc,'s','MarkerSize',12,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
    hold on;
    %errorbar(corrs,pc,upperCI-pc,lowerCI-pc); DON'T PLOT ERROR BARS
    % TOO FEW TRIALS TO PLOT ERROR BARS

    xlabel('SPLs'); ylabel('Percent Correct');
    %axis([-5 100 0 1]);
    title(['Summary Plot of ROC Analysis']);
    end
    
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
    %saveas(A.hFig,[A.fname(1:end-4) '_analysis.fig']);
    %saveas(A.hFig,[A.fname(1:end-4) '_analysis.png']);
    anal_fname=[A.PDRfile '_ANALYSIS.mat'];
    save([A.analysis_dir anal_fname],'A');
end
cd(current_dir);