% script and info to compare Env_Vstr for stimulus
% with Zscore Vstr of spikes
% implemented 10/11/04 for 964DD


% first get the stimulus summary produced earlier
load c:\kip_overflow\2sources\cues\abl\abl

% adjust the dir from the stimulus summary to reflect the center of the
% array
dir(1,:) = dir(1,:)-15;
dir(2,:) = dir(2,:)+15;

% read the spike files (change within the spet.. script to e.g. process
% only the abl test)
sept2804a

% set the spikes dir
dir_spikes = params.locations;

% find the indices into the two dirs (spikes and stim) to use matching locs
clear ind ind_s*
for i = 1:91
temp = find(dir_spikes(1,:) == dir(1,i) & dir_spikes(2,:) == dir(2,i));
if ~isempty(temp);
ind(i) = temp;
end
end
ind_stim = find(ind>0);
ind_spikes = ind(ind_stim);

% check out the dirs
plot_diam(Z_list(55,ind_spikes),dir(:,ind_stim),1);
plot_diam(Z_list(55,:),dir_spikes,1);

% run corrcoefs:
corrcoef(Z_list(55,ind_spikes),mean(Vstr_ABL_55(ind_stim,:),2))
corrcoef(Z_list(55,ind_spikes),mean(Vstr_dILD_55(ind_stim,:),2))
corrcoef(Z_list(55,ind_spikes),mean(Vstr_dIPD_55(ind_stim,:),2))
corrcoef(Z_list(75,ind_spikes),mean(Vstr_ABL_75(ind_stim,:),2))
corrcoef(Z_list(75,ind_spikes),mean(Vstr_dILD_75(ind_stim,:),2))
corrcoef(Z_list(75,ind_spikes),mean(Vstr_dIPD_75(ind_stim,:),2))
