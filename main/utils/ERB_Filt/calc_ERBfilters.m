function [fcoefs,Factor] = calc_ERBfilters;

% function to calculate ERB filters and normalization factors
% usually used in acoustics model(s)


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
Factor = Factor';