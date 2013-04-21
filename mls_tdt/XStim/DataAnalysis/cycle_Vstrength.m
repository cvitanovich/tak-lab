function [Vstr, Theta, Nspikes] = cycle_Vstrength(data, modPer)

%	function [Vstr, Theta, Nspikes] = cycle_Vstrength(data, modPer)
%	calculates Vector strength and best angle for
%	modulation periods in usec provided as argin
%	    - data should be a [nReps, nMaxSpikes] matrix of spiketimes in usec
%       - modPer in usec
%   - for use with McSpace and an AM that is locked to the epochs (e.g. 100
%   Hz atop 50 msec epochs)
%   - Vstr is calc'd for each modulation cycle (using reps to get
%   sufficient data)
%   assumes starting phase == 0

if nargin < 2
    disp('cycle_Vstrength requires two argins: data, modPer');
    return
end

ind = find(data < 0);
data(ind) = zeros(size(ind));

for icycle = 1:ceil(max1(data)/modPer)
    st = 1+(icycle-1)*modPer;
    en = icycle*modPer;
    ind = find(data>=st & data<en);
    if length(ind > 0)
        % calc phase for each spike
        phas = mod(data(ind),modPer) /modPer;   % in cycles
        phas = 2 * pi * phas;       % in radians
        % calc Vstr
        ycoord = mean(sin(phas(:)));
        xcoord = mean(cos(phas(:)));
        Vstr(icycle) = sqrt(ycoord^2 + xcoord^2);
        Nspikes(icycle) = length(ind);
        % calc best phase
        temp = atan2(ycoord,xcoord);
        if temp < 0 
            temp = temp + 2*pi;
        end
        Theta(icycle) = temp;
    else
        Vstr(icycle) = 0;
        Nspikes(icycle) = 0;
        Theta(icycle) = NaN;
    end
end