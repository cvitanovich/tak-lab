function [params,usedBank] = makeAMparams1(randBank,usedBank,irep)

% [params,usedBank] = makeAMparams1(randBank,usedBank,irep);
%
% randBank is [nReps x 5] a list of the seeds for each use of rand
% usedBank is [nReps x nTests] a list of which rows in randbank are used by
% each test (1 == used, 0 == available)
%
% makes parameters for Precedence Effect AM test, [#trials x 12]
%
% 1) seed for noise carrier
% 2) seed for decorrelating noise carrier
% 3) stimdur (ms)
% 4) delay (ms) (if delay==0, then uses decorrelated carriers)
% 5) ramp_dur (ms)
% 6) seed for AM envelope
% 7) gammatone cF (Hz)
% 8) envelope decorrelation index (0-1)
% 9) seed for decorrelating envelope#1
% 10) seed for decorrelating envelope#2
% 11) sound index for in RF   (1,2, or no Sound = 3)
% 12) sound index for in antiRF (1,2, or no Sound = 3)
%
%
% version 1: correcting use of usedBank 6/26/12
% will fail if irep > size(randBank,1)

global XStimParams


%randBank = round(rand(1,1000)*100000);
stimdur = XStimParams.curr_stimdur;
rampdur = 2.5;                  % should correspond to value in Engage_PE_env
delay = XStimParams.delay;
n_delay = length(delay);
cF = XStimParams.cF;
n_cF = length(cF);

% how many trials?
n_test = size(usedBank,2);
params = zeros(n_test,12);

itest = 1;

for iLoc=1:2
    for icF = 1:n_cF
        for idelay = 1:n_delay
            ind0 = setdiff(1:size(randBank,1),usedBank(:,itest)); % all available rows of randBank for this test
            ind1 = ind0(ceil(rand(1)*length(ind0)));    % which row to use for this test
            usedBank(irep,itest)=ind1;                   % mark it as used
            params(itest,1) = randBank(ind1,1);          % seed for noise carrier
            params(itest,2) = randBank(ind1,2);          % seed for decorrelating noise carrier
            params(itest,3) = stimdur;                  % stimdur (ms)
            params(itest,4) = delay(idelay);            % delay (ms) (if delay==0, then uses decorrelated carriers)
            params(itest,5) = rampdur;                  % ramp_dur (ms)
            params(itest,6) = randBank(ind1,3);          % seed for AM envelope
            params(itest,7) = cF(icF);                  % gammatone cF (Hz)
            params(itest,8) = 0;                        % envelope decorrelation index (0-1)
            params(itest,9) = randBank(ind1,4);          % seed for decorrelating envelope#1
            params(itest,10) = randBank(ind1,5);         % seed for decorrelating envelope#2
            switch iLoc
                case 1
                    params(itest,11) = 1;
                    params(itest,12) = 2;
                case 2
                    params(itest,11) = 2;
                    params(itest,12) = 1;
            end
            itest= itest+1;
        end
    end
end

% decorr env, decorr carrier
for icF = 1:n_cF
    ind0 = setdiff(1:size(randBank,1),usedBank(:,itest));         % all available rows of randBank for this test
    ind1 = ind0(ceil(rand(1)*length(ind0)));    % which row to use for this test
    usedBank(irep,itest)=ind1;                   % mark it as used
    params(itest,1) = randBank(ind1,1);          % seed for noise carrier
    params(itest,2) = randBank(ind1,2);          % seed for decorrelating noise carrier
    params(itest,3) = stimdur;                         % stimdur (ms)
    params(itest,4) = 0;                               % delay (ms) (if delay==0, then uses decorrelated carriers)
    params(itest,5) = rampdur;                         % ramp_dur (ms)
    params(itest,6) = randBank(ind1,3);                 % seed for AM envelope
    params(itest,7) = cF(icF);                         % gammatone cF (Hz)
    params(itest,8) = 1;                               % envelope decorrelation index (0-1)
    params(itest,9) = randBank(ind1,4);                 % seed for decorrelating envelope#1
    params(itest,10) = randBank(ind1,5);                % seed for decorrelating envelope#2
    params(itest,11) = 1;
    params(itest,12) = 2;
    itest= itest+1;
end

% identical env, decorr carrier
for icF = 1:n_cF
    ind0 = setdiff(1:size(randBank,1),usedBank(:,itest));         % all available rows of randBank for this test
    ind1 = ind0(ceil(rand(1)*length(ind0)));    % which row to use for this test
    usedBank(irep,itest)=ind1;                   % mark it as used
    params(itest,1) = randBank(ind1,1);          % seed for noise carrier
    params(itest,2) = randBank(ind1,2);          % seed for decorrelating noise carrier
    params(itest,3) = stimdur;                         % stimdur (ms)
    params(itest,4) = 0;                               % delay (ms) (if delay==0, then uses decorrelated carriers)
    params(itest,5) = rampdur;                         % ramp_dur (ms)
    params(itest,6) = randBank(ind1,3);                 % seed for AM envelope
    params(itest,7) = cF(icF);                         % gammatone cF (Hz)
    params(itest,8) = 0;                               % envelope decorrelation index (0-1)
    params(itest,9) = randBank(ind1,4);                 % seed for decorrelating envelope#1
    params(itest,10) = randBank(ind1,5);                % seed for decorrelating envelope#2
    params(itest,11) = 1;
    params(itest,12) = 2;
    itest= itest+1;
end

%%%%% these may have to have identical carriers?
% identical env, identical carrier one source in RF
for icF = 1:n_cF
    ind0 = setdiff(1:size(randBank,1),usedBank(:,itest));         % all available rows of randBank for this test
    ind1 = ind0(ceil(rand(1)*length(ind0)));    % which row to use for this test
    usedBank(irep,itest)=ind1;                   % mark it as used
    params(itest,1) = randBank(ind1,1);          % seed for noise carrier
    params(itest,2) = randBank(ind1,2);          % seed for decorrelating noise carrier
    params(itest,3) = stimdur;                         % stimdur (ms)
    params(itest,4) = 0;                               % delay (ms) (if delay==0, then uses decorrelated carriers)
    params(itest,5) = rampdur;                         % ramp_dur (ms)
    params(itest,6) = randBank(ind1,3);                 % seed for AM envelope
    params(itest,7) = cF(icF);                         % gammatone cF (Hz)
    params(itest,8) = 0;                               % envelope decorrelation index (0-1)
    params(itest,9) = randBank(ind1,4);                 % seed for decorrelating envelope#1
    params(itest,10) = randBank(ind1,5);                % seed for decorrelating envelope#2
    params(itest,11) = 1;
    params(itest,12) = 3;
    itest= itest+1;
end

% identical env, identical carrier one source in anti-RF
for icF = 1:n_cF
    ind0 = setdiff(1:size(randBank,1),usedBank(:,itest));         % all available rows of randBank for this test
    ind1 = ind0(ceil(rand(1)*length(ind0)));    % which row to use for this test
    usedBank(irep,itest)=ind1;                   % mark it as used
    params(itest,1) = randBank(ind1,1);          % seed for noise carrier
    params(itest,2) = randBank(ind1,2);          % seed for decorrelating noise carrier
    params(itest,3) = stimdur;                         % stimdur (ms)
    params(itest,4) = 0;                               % delay (ms) (if delay==0, then uses decorrelated carriers)
    params(itest,5) = rampdur;                         % ramp_dur (ms)
    params(itest,6) = randBank(ind1,3);                 % seed for AM envelope
    params(itest,7) = cF(icF);                         % gammatone cF (Hz)
    params(itest,8) = 0;                               % envelope decorrelation index (0-1)
    params(itest,9) = randBank(ind1,4);                 % seed for decorrelating envelope#1
    params(itest,10) = randBank(ind1,5);                % seed for decorrelating envelope#2
    params(itest,11) = 3;
    params(itest,12) = 1;
    itest= itest+1;
end