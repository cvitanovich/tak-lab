% script to analyse stimuli for 964 or 965 cells
% need to fill in many (%%%change%%%) variables below 

stim_dir = 'c:\kip_overflow\2sources\cues2\ITA2\';

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
Z_dILD_20 = zeros(nDir,n_cF);
Z_dILD_55 = zeros(nDir,n_cF);
Z_dILD_75 = zeros(nDir,n_cF);
Z_dIPD_20 = zeros(nDir,n_cF);
Z_dIPD_55 = zeros(nDir,n_cF);
Z_dIPD_75 = zeros(nDir,n_cF);
Z_ABL_20 = zeros(nDir,n_cF);
Z_ABL_55 = zeros(nDir,n_cF);
Z_ABL_75 = zeros(nDir,n_cF);

Vstr_dILD_20 = zeros(nDir,n_cF);
Vstr_dILD_55 = zeros(nDir,n_cF);
Vstr_dILD_75 = zeros(nDir,n_cF);
Vstr_dIPD_20 = zeros(nDir,n_cF);
Vstr_dIPD_55 = zeros(nDir,n_cF);
Vstr_dIPD_75 = zeros(nDir,n_cF);
Vstr_ABL_20 = zeros(nDir,n_cF);
Vstr_ABL_55 = zeros(nDir,n_cF);
Vstr_ABL_75 = zeros(nDir,n_cF);

for iDir = 1:nDir
    disp(['processing dir # ' num2str(iDir) ' of ' num2str(nDir)])
    
    % call function to calc the cues, envelopes, FTs, etc
    [ILD,ITD,ABL,cF,time,dIPD,dILD,Env_dIPD,Env_dILD,Env_ABL,FT_Env_dIPD,FT_Env_dILD,FT_Env_ABL] = ...
        plot_cues_from_mls_stims3([], 30, 0, best_IPD, ind_best_IPD, best_ILD,dir(1,iDir),dir(2,iDir),stim_dir);
    
    % get peaks of envelope FFTs
    [nFreq n] = size(FT_Env_ABL);
    freq = 0:500/(nFreq-1):500;
    ind20 = nearest_index(freq,20);
    ind55 = nearest_index(freq,55);
    ind75 = nearest_index(freq,75);
    
    for icF  = 1:n_cF
        Vstr_dILD_20(iDir,icF) = abs(FT_Env_dILD(ind20,icF));
        Vstr_dILD_55(iDir,icF) = max(abs(FT_Env_dILD([ind55 ind55+1],icF)));
        Vstr_dILD_75(iDir,icF) = abs(FT_Env_dILD(ind75,icF));
        temp = std(abs(FT_Env_dILD(2:nFreq,icF)));
        if temp >0
            Z_dILD_20(iDir,icF) = temp(ind20-1) / temp;
            Z_dILD_55(iDir,icF) = max(temp([ind55-1 ind55])) / temp;
            Z_dILD_75(iDir,icF) = temp(ind75-1) / temp;
        end
        
        Vstr_dIPD_20(iDir,icF) = abs(FT_Env_dIPD(ind20,icF));
        Vstr_dIPD_55(iDir,icF) = max(abs(FT_Env_dIPD([ind55 ind55+1],icF)));
        Vstr_dIPD_75(iDir,icF) = abs(FT_Env_dIPD(ind75,icF));
        temp = std(abs(FT_Env_dIPD(2:nFreq,icF)));
        if temp >0
            Z_dIPD_20(iDir,icF) = temp(ind20-1)/ temp;
            Z_dIPD_55(iDir,icF) = max(temp([ind55-1 ind55]))/ temp;
            Z_dIPD_75(iDir,icF) = temp(ind75-1)/ temp;
        end
        
        Vstr_ABL_20(iDir,icF) = abs(FT_Env_ABL(ind20,icF));
        Vstr_ABL_55(iDir,icF) = max(abs(FT_Env_ABL([ind55 ind55+1],icF)));
        Vstr_ABL_75(iDir,icF) = abs(FT_Env_ABL(ind75,icF));
        temp = std(abs(FT_Env_ABL(2:nFreq,icF)));
        if temp >0
            Z_ABL_20(iDir,icF) = temp(ind20-1)/ temp;
            Z_ABL_55(iDir,icF) = max(temp([ind55-1 ind55]))/ temp;
            Z_ABL_75(iDir,icF) = temp(ind75-1)/ temp;
        end
    end
end

%plot Vstrength with caxis normalized along rows
figure
m = min1([20*log10(mean(Vstr_ABL_20,2)) 20*log10(mean(Vstr_ABL_55,2)) 20*log10(mean(Vstr_ABL_75,2))])
M = max1([20*log10(mean(Vstr_ABL_20,2)) 20*log10(mean(Vstr_ABL_55,2)) 20*log10(mean(Vstr_ABL_75,2))])
% 20 Hz
h = axes('position',[.24 (.01 + (4)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(20*log10(mean(Vstr_ABL_20,2)), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (.01 + (4)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(20*log10(mean(Vstr_ABL_55,2)), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (.01 + (4)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(20*log10(mean(Vstr_ABL_75,2)), dir, 1);
caxis([m M])

m = min1([20*log10(mean(Vstr_dILD_20,2)) 20*log10(mean(Vstr_dILD_55,2)) 20*log10(mean(Vstr_dILD_75,2))])
M = max1([20*log10(mean(Vstr_dILD_20,2)) 20*log10(mean(Vstr_dILD_55,2)) 20*log10(mean(Vstr_dILD_75,2))])
% 20 Hz
h = axes('position',[.24 (.01 + (3)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(20*log10(mean(Vstr_dILD_20,2)), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (.01 + (3)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(20*log10(mean(Vstr_dILD_55,2)), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (.01 + (3)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(20*log10(mean(Vstr_dILD_75,2)), dir, 1);
caxis([m M])

m = min1([20*log10(mean(Vstr_dIPD_20,2)) 20*log10(mean(Vstr_dIPD_55,2)) 20*log10(mean(Vstr_dIPD_75,2))])
M = max1([20*log10(mean(Vstr_dIPD_20,2)) 20*log10(mean(Vstr_dIPD_55,2)) 20*log10(mean(Vstr_dIPD_75,2))])
% 20 Hz
h = axes('position',[.24 (.01 + (2)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(20*log10(mean(Vstr_dIPD_20,2)), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (.01 + (2)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(20*log10(mean(Vstr_dIPD_55,2)), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (.01 + (2)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(20*log10(mean(Vstr_dIPD_75,2)), dir, 1);
caxis([m M])

% plot Zscore Vstr with caxis normalized for all plots        
figure
temp = [mean(Z_ABL_20,2) mean(Z_ABL_55,2) mean(Z_ABL_75,2) ...
        mean(Z_dILD_20,2) mean(Z_dILD_55,2) mean(Z_dILD_75,2) ...
        mean(Z_dIPD_20,2) mean(Z_dIPD_55,2) mean(Z_dIPD_75,2)];
m = min1(temp)
M = max1(temp)
% 20 Hz
h = axes('position',[.24 (.01 + (4)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_ABL_20,2), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (.01 + (4)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_ABL_55,2), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (.01 + (4)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_ABL_75,2), dir, 1);
caxis([m M])

% 20 Hz
h = axes('position',[.24 (.01 + (3)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dILD_20,2), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (.01 + (3)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dILD_55,2), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (.01 + (3)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dILD_75,2), dir, 1);
caxis([m M])

% 20 Hz
h = axes('position',[.24 (.01 + (2)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dIPD_20,2), dir, 1);
caxis([m M])
% 55 Hz
h = axes('position',[.44 (.01 + (2)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dIPD_55,2), dir, 1);
caxis([m M])
% 75 Hz
h = axes('position',[.64 (.01 + (2)*.20) .21 .21]);
set(h,'XTickLabel',[])
set(h,'YTickLabel',[])
plot_diam_2axis(mean(Z_dIPD_75,2), dir, 1);
caxis([m M])
