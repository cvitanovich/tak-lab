function [] = anal_cycleVstr;

global H
global XStimParams
global FN

% anal_cycleVstr
% run analysis of Vstrength on cycle-by-cycle basis
% for data from McSpace using AMed stimuli locked to epochduration

set(H.mcSpace_analVstr,'value',0)

% check modtype and modfreq
if strfind(XStimParams.mod_type,'Tone') | strfind(XStimParams.mod_type,'Sq')
    if mod(XStimParams.epoch_duration(1), 1000/XStimParams.mod_freq(1))
        errordlg('analyse Vstr not allowed - modPer must be a factor of epochDuration')
        return;
    end
else
    errordlg('analyse Vstr not allowed - modType must be tone or sq_wave')
end

% query for filename
Prompt{1} = 'Path :';
Prompt{2} = 'Filename :';
Title = 'Enter data filename for vector analysis';
DefAns{1} = FN.data_path;
DefAns{2} = FN.data;
LineNo=1;
Answer = inputdlg(Prompt,Title,LineNo,DefAns);

% check for existance of this file
if ~exist([Answer{1} Answer{2} '.mat'],'file')
    errordlg('analyse Vstr not allowed - file does not exist')
    return;
end    

% load DATA and change data to usecs
eval(['load ' Answer{1} Answer{2} ' DATA params param3'])

data = DATA(2:end,DATA(1,1)+2:end);
data = data - params.silence_lead;
data = 1000*(data - 3.1)*1.0016;
ind = find(data<0);
data(ind) = zeros(size(ind));

% run analysis and plot
[Vstr, Theta, Nspikes] = cycle_Vstrength(data, 1000000/XStimParams.mod_freq(1));

figure; hold on

cycles = 1:1000000/params.mod_freq(1):params.curr_stimdur*1000;
Ncycles = length(cycles);
if length(Vstr)<Ncycles
    Vstr = [Vstr zeros(1,Ncycles-length(Vstr))];
end
plot(cycles, Vstr(1:Ncycles),'b.');

cyc_per_epoch = params.epoch_duration(1) / (1000/params.mod_freq(1));
levels = zeros(size(cycles));
for icycle = 1:cyc_per_epoch
    levels(icycle:cyc_per_epoch:end) = param3(1,:)/max1(param3(1,:));
end
plot(cycles,levels,'r');