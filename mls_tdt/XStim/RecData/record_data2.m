function [] = record_data2(params,datamatrix)
%Spike Data Recording Routine version2
% useful for spatial data

% separate spikes by trial
% find maxspikes
i = 0; j = 1;
while i < size(datamatrix,1)
   Nspikes(j) = datamatrix(i+1,1);
   i = i + Nspikes(j)+1;
   j = j+1;
end

% make DATA array (el az rep Nspikes spiketimes...)
maxSpikes = max1(Nspikes);
DATA = zeros(j,maxSpikes+4);

i = 0; j = 1;
while i < size(datamatrix,1)
   DATA(j,1:Nspikes(j) + 4) = [datamatrix(i+1,3:4) datamatrix(i+1,2) ...
         Nspikes(j) datamatrix(i+2:i+Nspikes(j)+1)];
   i = i + Nspikes(j)+1;
	j = j+1;   
end


datafile = [params.datadir params.datafile];
eval(['save ' datafile ' params DATA']);

return;


% DATA:
%	el	az	repnum	Nspikes	spiketimes...