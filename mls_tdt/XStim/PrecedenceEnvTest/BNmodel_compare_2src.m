function [SP1, SP2, leadP, lagP] = BNmodel_compare_2src(Env1, Env2)

% [SP1, SP2, leadP, lagP] = BNmodel_compare_2src(Env1, Env2);
%
% Brian's model for owl stimuli
% inputs are linear envelopes
% outputs are spike probabilities as function of time and summed

Fs = 30000;

% get amplitude difference
Env1(Env1==0) = 0.00000000001;
Env1_dB = 20*log10(Env1);
Env2(Env2==0) = 0.00000000001;
Env2_dB = 20*log10(Env2);
Env1_dBDiff = Env1_dB - Env2_dB;
Env2_dBDiff = Env2_dB - Env1_dB;%-Env1_dBDiff;
Env1_dBDiff_wt = 3.0207 ./(1 + exp(-(Env1_dBDiff(:) - 3.2497) ./ 2.3835));
Env2_dBDiff_wt = 3.0207 ./(1 + exp(-(Env2_dBDiff(:) - 3.2497) ./ 2.3835));
%Env1_dBDiff = max(-24, Env1_dBDiff); % limit range
%Env2_dBDiff = min(24, Env2_dBDiff); % limit range

% derivative
% use average envelope to measure d/dt -- as in Nelson & Takahashi 2010
% calculate average envlope
Env1_dt = diff((Env1 + Env2)/2) * (Fs/1000); % units are per ms... from AM paper
Env2_dt = diff((Env1 + Env2)/2) * (Fs/1000);

if 0
    Env1_dt_wt = 2.2082 ./(1 + exp(-(Env1_dt(:) - 0.014001) ./ 0.065264));
    Env2_dt_wt = 2.2082 ./(1 + exp(-(Env2_dt(:) - 0.014001) ./ 0.065264));
else                % shifts the weighting function to the right
    Env1_dt_wt = 2.2082 ./(1 + exp(-(Env1_dt(:) - 0.2) ./ 0.065264));
    Env2_dt_wt = 2.2082 ./(1 + exp(-(Env2_dt(:) - 0.2) ./ 0.065264));
end

SP1 = (Env1_dt_wt) .* (Env1_dBDiff_wt(2:end));
SP2 = (Env2_dt_wt) .* (Env2_dBDiff_wt(2:end));

if 1 % entire stimulus
    leadP = sum(SP1) / (Fs/1000) / 34; % 34=duration of stimuli used to estimate sensitivity functions?
    lagP = sum(SP2) / (Fs/1000) / 34;
else % skip onset and offset
    firstpnt = round(0.003*Fs);
    lastpnt = round(length(SP1)-(0.003*Fs));
    leadP = sum(SP1(firstpnt:lastpnt)) / (Fs/1000) / 34;
    lagP = sum(SP2(firstpnt:lastpnt)) / (Fs/1000) / 34;
end


