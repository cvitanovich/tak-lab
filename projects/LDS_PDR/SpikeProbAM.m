function [leadP, lagP] = SpikeProbAM(Env1, Env2,Fs)
% Adapted from B.S. Nelson's AM study
% leadP = lead spike probability (cummulative)
% lagP = lag spike probability (cummulative)

% Env1 = lead envelope
% Env2 = lag envelope

SR = Fs; % may need to adjust derivative function if SR is changed
% set separately in function filteredstimuliAM !!!

% get amplitude difference
Env1(find(Env1==0)) = 0.00000000001;
Env1_dB = 20*log10(Env1);
Env2(find(Env2==0)) = 0.00000000001;
Env2_dB = 20*log10(Env2);
Env1_dBDiff = Env1_dB - Env2_dB;
Env2_dBDiff = Env2_dB - Env1_dB;%-Env1_dBDiff;
Env1_dBDiff_wt = 3.0207 ./(1 + exp(-(Env1_dBDiff(:) - 3.2497) ./ 2.3835));
Env2_dBDiff_wt = 3.0207 ./(1 + exp(-(Env2_dBDiff(:) - 3.2497) ./ 2.3835));

%Env1_dBDiff = max(-24, Env1_dBDiff); % limit range
%Env2_dBDiff = min(24, Env2_dBDiff); % limit range

% derivative
if 1 % use average envelope to measure d/dt -- as in Nelson & Takahashi 2010
    % calculate average envlope
    Env1_dt = diff((Env1 + Env2)/2) * (SR/1000); % units are per ms... from AM paper
    Env2_dt = diff((Env1 + Env2)/2) * (SR/1000);
   
    Env1_dt_wt = 2.2082 ./(1 + exp(-(Env1_dt(:) - 0.014001) ./ 0.065264));
    Env2_dt_wt = 2.2082 ./(1 + exp(-(Env2_dt(:) - 0.014001) ./ 0.065264));
   
    spikeProb_1 = (Env1_dt_wt) .* (Env1_dBDiff_wt(2:end));
    spikeProb_2 = (Env2_dt_wt) .* (Env2_dBDiff_wt(2:end));
   
else % calulate d/dt separately for each sound -- as in human model
    % This doesn't work (here) unless we get signals spatialized in
    % each ear. The derivavite would then be calulated from the
    % signals in each ear.
    % Also cannot use functions from Nelson & Takahashi 2010
    Env1_dt = diff((Env1)/2) * (SR/1000); % units are per ms... from AM paper
    Env2_dt = diff((Env2)/2) * (SR/1000);
    Env1_dt_wt = 2.2082 ./(1 + exp(-(Env1_dt(:) - 0.014001) ./ 0.065264));
    Env2_dt_wt = 2.2082 ./(1 + exp(-(Env2_dt(:) - 0.014001) ./ 0.065264));
    Env_dt_wt = Env1_dt_wt .* Env2_dt_wt; % multiply weights
    Env1_dt_wt = Env_dt_wt;
    Env2_dt_wt = Env_dt_wt;
    spikeProb_1 = (Env_dt_wt) .* (Env1_dBDiff_wt(2:end));
    spikeProb_2 = (Env_dt_wt) .* (Env2_dBDiff_wt(2:end));
end

if 1 % entire stimulus
    leadP = sum(spikeProb_1) / (SR/1000) / 34; % 34=duration of stimuli used to estimate sensitivity functions?
    lagP = sum(spikeProb_2) / (SR/1000) / 34;
else % skip onset and offset
    firstpnt = round(0.003*SR);
    lastpnt = round(length(spikeProb_1)-(0.003*SR));
    leadP = sum(spikeProb_1(firstpnt:lastpnt)) / (SR/1000) / 34;
    lagP = sum(spikeProb_2(firstpnt:lastpnt)) / (SR/1000) / 34;
end