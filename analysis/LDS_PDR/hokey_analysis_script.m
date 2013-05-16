% LOAD LIST OF FILES
anal_dir='/Users/cvitanovich/Documents/MATLAB/analysis_data';
data_dir='/Users/cvitanovich/Desktop/alex_data/';
list_dir='/Users/cvitanovich/Documents/MATLAB/tak-lab/analysis/LDS_PDR';
ZTRIALS=[];
ZHAB=[];
cd(list_dir)
fnames=textread('file_list.txt','%s',-1);
cd(anal_dir)
nFiles=size(fnames,1);
hWait=waitbar(0,'grabbing z-scores...');
for g=1:nFiles
    % LOAD ANALYSIS FILE
    clear A PDR
    F=[fnames{g} '_ANALYSIS.mat'];
    load(F);
    % LOAD HEADER FILE
    H=[data_dir fnames{g} '.mat'];
    load(H);
    delay=PDR.SOUNDS_carrier_delay;
    % For each loc,...
    nTypes=length(A.trials.types);
    for h=1:nTypes
        % azimuthal shift ("DELTA")
        DELTA = str2num(A.trials.test{h}.loc)-PDR.LAG_hab_pos;
        % For each trial, grab z-scores and add to matrix:
        Z=A.trials.test{h}.zscores;
        nTrials=length(Z);
        for j=1:nTrials
            % store trial in combined matrix!
            ZTRIALS=[ZTRIALS; DELTA Z(j) delay];
        end
    end
    % grab z-scores for hab trials!
    Zes=A.trials.hab.zscores';
    Des=zeros(length(Zes),1);
    Delays=delay*ones(length(Zes),1);
    ZHAB=[ZHAB; Des Zes Delays];
    waitbar(g/nFiles,hWait);
end
close(hWait)

% plot z-scores!!!!
figure; hold on;
nTrials=size(ZTRIALS,1);
hWait=waitbar(0,'plotting...');
for k=1:nTrials
    if(ZTRIALS(k,3)==-3)
        col='r';
    else
        col='g';
    end
    scatter(ZTRIALS(k,1),ZTRIALS(k,2),10,col);
    waitbar(k/nTrials,hWait)
end
close(hWait);

% standard separation:
std_sep=@(m1,m2,s1,s2) (m1-m2)/sqrt(s1*s2);

figure; hold on;
Des=[];
DELTAS=unique(ZTRIALS(:,1));
for l=1:length(DELTAS)
    
    [row0 col]=find(ZTRIALS(:,1)==DELTAS(l));
    [row1 col]=find(ZTRIALS(:,3)==-3);
    idx=intersect(row0,row1);
    zlead=ZTRIALS(idx,2); % z-scores for delay = -3ms, this delta!
    [row2 col]=find(ZTRIALS(:,3)==3);
    idx=intersect(row0,row2);
    zlag=ZTRIALS(idx,2); % z-scores for delay = +3ms, this delta!
    %hab trials:
    [row col]=find(ZHAB(:,3)==-3);
    zhab_lead=ZHAB(row,2);
    [row col]=find(ZHAB(:,3)==3);
    zhab_lag=ZHAB(row,2);
    
    % standard sep. (lead, this delta)
    mlead=mean(zlead);
    slead=std(zlead);
    mhab_lead=mean(zhab_lead);
    shab_lead=std(zhab_lead);
    Dlead=std_sep(mlead,mhab_lead,slead,shab_lead);
    
    % standard separation (lag, this delta)
    mlag=mean(zlag);
    slag=std(zlag);
    mhab_lag=mean(zhab_lag);
    shab_lag=std(zhab_lag);
    Dlag=std_sep(mlag,mhab_lag,slag,shab_lag);
    
    % plotting
    scatter(DELTAS(l),Dlead,20,'r');
    scatter(DELTAS(l),Dlag,20,'g');
    Des=[Des; DELTAS(l) Dlead Dlag];
end

    
    