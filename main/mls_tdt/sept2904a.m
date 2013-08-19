figure
pcolor((1:nWin),cF,abs(dILD_p40n50)')
    shading flat
    title('dILD p40 n50')
hold on
plot((1:100),mean(ABL_p40n50,2)*200,'k','linewidth',3)


figure; hold on
plot(ABL_p40n50 * weights, 'b','linewidth',2)
plot(ABL_p10n20 * weights, 'r','linewidth',2)
plot((ABL_p40n50 * weights + ABL_p10n20 * weights)/2, 'k','linewidth',2)

plot((abs(dILD_p40n50) * weights + 30)*2.5, 'c','linewidth',2)
plot((abs(dILD_p10n20) * weights + 30)*2.5, 'm','linewidth',2)




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
Factormat = repmat(Factor',1,win);
clear S tempL M
Factormat = repmat(Factor',1,255);
tempL = ERBFilterBankB(TF1(308,:), fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
tempR = ERBFilterBankB(TF2(308,:), fcoefs) .* Factormat;
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
        
        
        
        