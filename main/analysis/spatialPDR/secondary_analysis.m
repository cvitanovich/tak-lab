% secondary analysis
h=figure; hold on;

% list of expt sessions
session_list={''};
pname='/home/andrew/analysis_data/';

%% MAIN LOOP
HabZes=[]; SPLs=[]; Zes=[];

for k=1:length(session_list)
    % load each session
    load([pname session_list{k} '_ANALYSIS.mat']);
    % grab hab z-scores
    HabZes=[HabZes A.trials.hab.zscores];
    % grab z-scores and spls
    for tt=1:length(A.trials.types)
        tmp=A.trials.test{tt}.zscores;
        SPLs=[SPLs str2num(A.trials.test{tt}.spl)*ones(1,length(tmp))];
        Zes=[Zes tmp];
    end
end

SPL_list=unique(SPLs); % list of unique SPLs tested

%% PARSE and analyze together

% ROC analysis by spl
warning off stats:perfcurve:SubSampleWithMissingClasses
dim=ceil(sqrt(length(SPL_list)));
for s=1:length(SPL_list)
    level=SPL_list(s);
    idx=find(SPLs==level);
    ztest=Zes(idx);
    zhabit = randsample(HabZes,100);
    scores = [zhabit ztest]; % combined areas
    labels(1:length(zhabit)) = 0;
    labels(length(zhabit)+1:length(scores)) = 1;
    [X,Y,T,AUC] = perfcurve(labels,scores,1,'NBoot',10);
    figure(h); hold on;
    hSub=subplot(dim,dim,s);
    %X2=[X(:,1)'; X(:,1)'; X(:,1)'];
    %Y2=[Y(:,2)' - Y(:,1)'; Y(:,1)'; Y(:,3)' - Y(:,1)'];
    h=area(Y(:,1)'); set(h,'FaceColor',[s/length(SPL_list) 0 0]);
    %colormap summer
    %set(h,'FaceColor',[i/length(A.trials.types) 0 0]);
    %set(h,'Marker','none','Color',colr,'LineWidth',3); %'LineStyle',style_type);
    % plot errors
    %errorbar(X(:,1),Y(:,1),Y(:,2)-Y(:,1),Y(:,3)-Y(:,1));
    axis([1 length(Y(:,1)') 0 1]);
    len=length(Y(:,1)');
    set(hSub,'XTick',1:(len/2):len,'XTickLabel',{'0','0.5','1'});
    xlabel('FPR'); ylabel('TPR');
    pc(i) = AUC(1);
    lowerCI(i) = AUC(2);
    upperCI(i) = AUC(3);
    title(['SPL = ' num2str(level) ' dB, PC=' num2str(AUC(1))]);
    drawnow;
end