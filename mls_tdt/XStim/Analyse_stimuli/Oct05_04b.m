% script to analyse stimuli for 964 or 965 cells
% need to fill in many (%%%change%%%) variables below 

cue = 'ITA2';
cellID = '965FA';
stim_dir = ['e:\kip\temp\' cue '\'];
sflag = 0;           %save flag
% initial params:
Fs = 30000;
cF = round(1000*exp(([12:40]/12)*log(2)))';
n_cF = length(cF);

% calc ERBfilters
[fcoefs,Factor] = calc_ERBfilters;
Factormat = repmat(Factor',1,255);

%%%%%%%%%%%%%%%%%%% find best ILD
dir = 0;
%%%change%%%
load -mat e:\kip\hrtfdata\965af.ila.mat   
best_EL = 5;
best_AZ = 5;

ind = find(dir(1,:) == best_EL & dir(2,:) == best_AZ);
tempL = ERBFilterBankB(TF1(ind,:), fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
tempR = ERBFilterBankB(TF2(ind,:), fcoefs) .* Factormat;
best_ILD = calclevel_time(tempL,tempR, cF);
clear TF* 

%%%%%%%%%%%%%%%%% find best IPD and ind_best_IPD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir = 0;
%%%change%%%
load -mat e:\kip\hrtfdata\965af.ita.mat   
best_EL = 5;
best_AZ = 5;

ind = find(dir(1,:) == best_EL & dir(2,:) == best_AZ);
[IPD, ind_best_IPD] = calc_bestIPD(TF1(ind,:),TF2(ind,:));
% get the IPD max
for icF = 1:n_cF
    best_IPD(icF) = IPD(icF,ind_best_IPD(icF));
end
clear IPD TF HDR dir

% start loop through dir here
dir = 0;
%%%change%%%
eval(['load ' stim_dir 'dir']);

dir = dir';
dir(1,:) = dir(1,:)+15;
dir(2,:) = dir(2,:)-15;
[n nDir] = size(dir);

% initialize arrays
nFreq = 256;
Z_dILD = zeros(nDir,n_cF,nFreq);
Z_dIPD = zeros(nDir,n_cF,nFreq);
Z_ABL = zeros(nDir,n_cF,nFreq);

Vstr_dILD = zeros(nDir,n_cF,nFreq);
Vstr_dIPD = zeros(nDir,n_cF,nFreq);
Vstr_ABL = zeros(nDir,n_cF,nFreq);

M_dILD = zeros(nDir,n_cF);
M_dIPD = zeros(nDir,n_cF);
M_ABL = zeros(nDir,n_cF);

for iDir = 1:nDir
    disp(['processing dir # ' num2str(iDir) ' of ' num2str(nDir)])
    
    % call function to calc the cues, envelopes, FTs, etc
    [ILD,ITD,ABL,cF,time,dIPD,dILD,Env_dIPD,Env_dILD,Env_ABL,FT_Env_dIPD,FT_Env_dILD,FT_Env_ABL] = ...
        plot_cues_from_mls_stims3([], 30, 0, best_IPD, ind_best_IPD, best_ILD,dir(1,iDir),dir(2,iDir),stim_dir);
    
    M_dILD(iDir,:) = mean(dILD(3:end,:),1);
    M_dIPD(iDir,:) = mean(dIPD(3:end,:),1);
    M_ABL(iDir,:) = mean(ABL(3:end,:),1);
    
    for icF  = 1:n_cF
        Vstr_dILD(iDir,icF,:) = abs(FT_Env_dILD(:,icF))';
        Vstr_dIPD(iDir,icF,:) = abs(FT_Env_dIPD(:,icF))';
        Vstr_ABL(iDir,icF,:) = abs(FT_Env_ABL(:,icF))';
    end
    
    for icF = 1:n_cF
        %disp(['processing cF # ' num2str(icF) ' of ' num2str(n_cF)])
        temp = Vstr_dILD(iDir,icF,:);
        ind = find(temp>0);
        temp = temp(ind);
        temp = std(temp(:));
        Z_dILD(iDir,icF,:) = Vstr_dILD(iDir,icF,:)/temp;
        
        temp = Vstr_dIPD(iDir,icF,:);
        ind = find(temp>0);
        temp = temp(ind);
        temp = std(temp(:));
        Z_dIPD(iDir,icF,:) = Vstr_dIPD(iDir,icF,:)/temp;
        
        temp = Vstr_ABL(iDir,icF,:);
        ind = find(temp>0);
        temp = temp(ind);
        temp = std(temp(:));
        Z_ABL(iDir,icF,:) = Vstr_ABL(iDir,icF,:)/temp;
    end
end
% get peaks of envelope FFTs
freq = 0:500/(nFreq-1):500;
ind20 = nearest_index(freq,20);
ind55 = nearest_index(freq,55);
ind75 = nearest_index(freq,75);


% plot
figure('position',[10 50 1000 650]);
text(.9, 1.05,'ABL2','fontsize',12)
text(.9, .8,'ILA2','fontsize',12)
text(.9, .55,'ITA2','fontsize',12)
text(.9, .25,'ABL2','fontsize',12)
text(.9, .00,'ILA2','fontsize',12)
text(.25, .35,'20 Hz                                    55 Hz                                        75 Hz','fontsize',12)
text(-.1, 1.05, cue,'fontsize',18)
text(-.1, 1.0, cellID,'fontsize',18)

text(-.1, .90, 'Zscore FT-Vstrength','fontsize',18)
text(-.1, .85, 'from stimuli','fontsize',18)
text(-.1, .65, 'caxis normalized','fontsize',18)
text(-.1, .60, 'throughout','fontsize',18)
text(-.1, .37, '- - - - - - - - - - - - - ','fontsize',18)
text(-.1, .20, 'caxis normalized','fontsize',18)
text(-.1, .15, 'by row','fontsize',18)
axis off

temp = [mean(Z_ABL(:,:,ind20),2) mean(Z_ABL(:,:,ind55),2) mean(Z_ABL(:,:,ind75),2) ...
        mean(Z_dILD(:,:,ind20),2) mean(Z_dILD(:,:,ind55),2) mean(Z_dILD(:,:,ind75),2) ...
        mean(Z_dIPD(:,:,ind20),2) mean(Z_dIPD(:,:,ind55),2) mean(Z_dIPD(:,:,ind75),2)];
m = min1(temp);
M = max1(temp);
%% ABL
% 20 Hz
h = axes('position',[.24 (.01 + (4)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_ABL(:,:,ind20),2), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (.01 + (4)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_ABL(:,:,ind55),2), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (.01 + (4)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_ABL(:,:,ind75),2), dir, 1);
caxis([m M])

%dILD
% 20 Hz
h = axes('position',[.24 (.01 + (3)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dILD(:,:,ind20),2), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (.01 + (3)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dILD(:,:,ind55),2), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (.01 + (3)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dILD(:,:,ind75),2), dir, 1);
caxis([m M])

%dIPD
% 20 Hz
h = axes('position',[.24 (.01 + (2)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dIPD(:,:,ind20),2), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (.01 + (2)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dIPD(:,:,ind55),2), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (.01 + (2)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dIPD(:,:,ind75),2), dir, 1);
caxis([m M])

% repeat with row-wise caxis scaling     
%% ABL
temp = [mean(Z_ABL(:,:,ind20),2) mean(Z_ABL(:,:,ind55),2) mean(Z_ABL(:,:,ind75),2)];
m = min1(temp);
M = max1(temp);
% 20 Hz
h = axes('position',[.24 (.01 + .15) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_ABL(:,:,ind20),2), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (.01 + .15) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_ABL(:,:,ind55),2), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (.01 + .15) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_ABL(:,:,ind75),2), dir, 1);
caxis([m M])

%dILD
temp = [mean(Z_dILD(:,:,ind20),2) mean(Z_dILD(:,:,ind55),2) mean(Z_dILD(:,:,ind75),2)];
m = min1(temp);
M = max1(temp);
h = axes('position',[.24 (-.05) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dILD(:,:,ind20),2), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (-.05) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dILD(:,:,ind55),2), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (-.05) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dILD(:,:,ind75),2), dir, 1);
caxis([m M])

%%%% save
if sflag
eval(['save e:\kip\temp\' cue '\' cue]);
end