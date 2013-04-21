% rasterize_DATA from read_data3_minus_spont

Nparams = DATA(1,1);
Nreps = params.numreps;
maxspikes = max1(DATA(:,Nparams+1));
frame = [0 params.curr_stimdur];
spiketimes = [];

figure; hold on;
for trial = 2:size(DATA,1)
       rep = DATA(trial,Nparams);
       ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
       ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(1)));
       ind2 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(2)));
       ind = setdiff(ind1,ind2);
       spiketimes = [spiketimes DATA(trial,ind + Nparams+1)];
       plot(DATA(trial,ind + Nparams+1),trial * ones(size(ind)),'.');
end         % end of trial

N = histc(spiketimes,[100:2:200]);
figure
bar([100:2:200],N,'histc');