% script to get best ILD and best IPD for 964 or 965 cells
% need to fill in many (%%%) variables below


% calc ERBfilters
Fs = 30000;
% calc cF (center frequencies for filterbank) on 1/12th octave scale
cF = round(1000*exp(([12:40]/12)*log(2)))';
n_cF = length(cF);
fcoefs = Make_ERBFiltA(Fs,cF);

S = zeros(1024*16,1);
S(1024*8) = 1;
temp = ERBFilterBankB(S, fcoefs);
for icF = 1:n_cF
    M(icF) = mom(temp(icF,:),2);
end
Factor = max1(M)./M;
Factormat = repmat(Factor',1,255);
clear S tempL M



%%%%%%%%%%%%%%%%%%% find best ILD
dir = 0;
load -mat e:\kip\hrtfdata\965af.ila.mat   
best_EL = -20;
best_AZ = 0;

ind = find(dir(1,:) == best_EL & dir(2,:) == best_AZ);
tempL = ERBFilterBankB(TF1(ind,:), fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
tempR = ERBFilterBankB(TF2(ind,:), fcoefs) .* Factormat;
best_ILD = calclevel_time(tempL,tempR, cF);

clear TF* 


%%%%%%%%%%%%%%%%% find best IPD and ind_best_IPD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir = 0;
load -mat e:\kip\hrtfdata\965af.ita.mat   
best_EL = -20;
best_AZ = 0;

ind = find(dir(1,:) == best_EL & dir(2,:) == best_AZ);
tempL = ERBFilterBankB(TF1(ind,:), fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
tempR = ERBFilterBankB(TF2(ind,:), fcoefs) .* Factormat;
[best_ITD, time, IPD] = calcitd(tempL,tempR, cF, Fs, ones(size(cF)));
temp = IPD(:,201:350);          % just the central portion (minimizes problems of phase ambiguity)
for icF = 1:n_cF
    zero = find(temp(icF,:) == max1(temp(icF,:))) + 200;
    plus = nearest_index(time,time(zero) + 1/cF(icF));
    minus = nearest_index(time,time(zero) - 1/cF(icF));
    plus2 = nearest_index(time,time(zero) + 2/cF(icF));
    minus2 = nearest_index(time,time(zero) - 2/cF(icF));
    if icF == 1
        ind_best_IPD(icF) = zero;
    else
        d0 = [zero minus plus minus2 plus2];
        d = abs(d0 - ind_best_IPD(icF-1));
        d1 = find(d == min1(d));
        ind_best_IPD(icF) = d0(d1);
    end
end
% get the IPD max
for icF = 1:n_cF
    best_IPD(icF) = IPD(icF,ind_best_IPD(icF));
end

clear ITD IPD time TF* d* m* ind temp* zero HDR F*


% start loop through dir here
% dir comprises three diagonals
%dir = [0 0; 0 -20; 0 -40; -10 10; 10 -10; -10 -10; 10 -30; -10 -30; 10 -50; -15 15; 15 -15; -15 -25; ...
%        15 -35; -15 -5; 15 -55; -20 0; 20 -20; -20 -20; 20 -40; 20 -60; -25 -15; 25 -25; 25 -45; -25 5; 25 -65; -30 -10; ...
%        30 -30; 30 -50; 35 -35; -35 -5; 35 -55; 40 -40; 45 -45; -5 -15; 5 -25; -5 -35; 5 -45; -5 5; 5 -5]';
% three strips: sum to -40 or -20 or 0
%ind = find(dir(1,:) + dir(2,:) == -20);
%dir = dir(:,ind);


dir = [35   -25
    30   -20
    25   -15
    20   -10
    15    -5
    10     0
     5     5
     0    10
    -5    15
   -10    20
   -15    25
   -20    30
   -25    35
    30   -30
    25   -25
    20   -20
    15   -15
    10   -10
     5    -5
     0     0
    -5     5
   -10    10
   -15    15
   -20    20
   -25    25
   -30    30
    25   -35
    20   -30
    15   -25
    10   -20
     5   -15
     0   -10
    -5    -5
   -10     0
   -15     5
   -20    10
   -25    15
   -30    20
   -35    25
    20   -40
    15   -35
    10   -30
     5   -25
     0   -20
    -5   -15
   -10   -10
   -15    -5
   -20     0
   -25     5
   -30    10
   -35    15
   -40    20
    15   -45
    10   -40
     5   -35
     0   -30
    -5   -25
   -10   -20
   -15   -15
   -20   -10
   -25    -5
   -30     0
   -35     5
   -40    10
   -45    15
    10   -50
     5   -45
     0   -40
    -5   -35
   -10   -30
   -15   -25
   -20   -20
   -25   -15
   -30   -10
   -35    -5
   -40     0
   -45     5
   -50    10
     5   -55
     0   -50
    -5   -45
   -10   -40
   -15   -35
   -20   -30
   -25   -25
   -30   -20
   -35   -15
   -40   -10
   -45    -5
   -50     0
   -55     5]';


load e:\kip\temp\ABA2\dir
dir = dir';
dir(1,:) = dir(1,:)+15;
dir(2,:) = dir(2,:)-15;


[n nDir] = size(dir);

for iDir = 1:nDir
    disp(['processing dir # ' num2str(iDir) ' of ' num2str(nDir)])
    
    % call function to calc the cues, envelopes, FTs, etc
    [ILD,ITD,ABL,cF,time,dIPD,dILD,Env_dIPD,Env_dILD,Env_ABL,FT_Env_dIPD,FT_Env_dILD,FT_Env_ABL] = ...
        plot_cues_from_mls_stims3([], 30, 0, best_IPD, ind_best_IPD, best_ILD,dir(1,iDir),dir(2,iDir));
    
    % get peaks of envelope FFTs
    [m n] = size(FT_Env_ABL);
    freq = 0:500/(m-1):500;
    ind20 = nearest_index(freq,20);
    ind55 = nearest_index(freq,55);         % actual freq is very close to mean(ind55 ind55+1)
    ind75 = nearest_index(freq,75);
    
    for icF  = 1:n_cF
        Vstr_dILD_20(iDir,icF) = abs(FT_Env_dILD(ind20,icF));
        Vstr_dILD_55(iDir,icF) = max1(abs(FT_Env_dILD([ind55 ind55+1],icF)));
        Vstr_dILD_75(iDir,icF) = abs(FT_Env_dILD(ind75,icF));
        
        Vstr_dIPD_20(iDir,icF) = abs(FT_Env_dIPD(ind20,icF));
        Vstr_dIPD_55(iDir,icF) = max1(abs(FT_Env_dIPD([ind55 ind55+1],icF)));
        Vstr_dIPD_75(iDir,icF) = abs(FT_Env_dIPD(ind75,icF));
        
        Vstr_ABL_20(iDir,icF) = abs(FT_Env_ABL(ind20,icF));
        Vstr_ABL_55(iDir,icF) = max1(abs(FT_Env_ABL([ind55 ind55+1],icF)));
        Vstr_ABL_75(iDir,icF) = abs(FT_Env_ABL(ind75,icF));
    end
end

%%%%%% add offsets to dir for plotting
%dir(1,:) = dir(1,:)+15;
%dir(2,:) = dir(2,:)-15;

% plot average (across freqs) Vstr for each cue and mod freq
plot_diam(20*log10(mean(Vstr_dILD_20,2)),dir);
title('dILD 20 Hz')
axis([-90 90 -90 90])
plot_diam(20*log10(mean(Vstr_dILD_55,2)),dir);
title('dILD 55 Hz')
axis([-90 90 -90 90])
plot_diam(20*log10(mean(Vstr_dILD_75,2)),dir);
title('dILD 75 Hz')
axis([-90 90 -90 90])

plot_diam(20*log10(mean(Vstr_dIPD_20,2)),dir);
title('dIPD 20 Hz')
axis([-90 90 -90 90])
plot_diam(20*log10(mean(Vstr_dIPD_55,2)),dir);
title('dIPD 55 Hz')
axis([-90 90 -90 90])
plot_diam(20*log10(mean(Vstr_dIPD_75,2)),dir);
title('dIPD 75 Hz')
axis([-90 90 -90 90])

plot_diam(20*log10(mean(Vstr_ABL_20,2)),dir);
title('ABL 20 Hz')
axis([-90 90 -90 90])
plot_diam(20*log10(mean(Vstr_ABL_55,2)),dir);
title('ABL 55 Hz')
axis([-90 90 -90 90])
plot_diam(20*log10(mean(Vstr_ABL_75,2)),dir);
title('ABL 75 Hz')
axis([-90 90 -90 90])
