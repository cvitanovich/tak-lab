function [DATA] = record_data3(params,datamatrix,param3,param4,param5,saveFlag)
%Spike Data Recording Routine version3		4/01/02
% datamatrix should come in as a list of spiketimes, 
% each associated with however many params are necessary
% to fuly describe it:
% Nspikes this rep	spiketime1	param1	param2	param3...
% Nspikes this rep	spiketime2	param1	param2	param3...
% Nspikes this rep	spiketime3	param1	param2	param3...
%% data will be saved as:
%	param1	param2	param3 ...	repnum	Nspikes	spiketimes...
% first line is a dummy with the number of params (including repnum) written to param1
% added param3 as argin to keep track of randomization for each rep
%
% added saveFlag to allow computation of DATA without saving

global FN

% find Ntrials and maxspikes
i = 1; Ntrials = 0; Nspikes = [];
while i < size(datamatrix,1)
   Nspikes(Ntrials+1) = datamatrix(i,1);
   i = i + Nspikes(Ntrials+1);
   Ntrials = Ntrials+1;
end
maxSpikes = max1(Nspikes);

Nparams = size(datamatrix,2)-2;
DATA = zeros(Ntrials+1,maxSpikes+Nparams+1);
DATA(1,1) = Nparams;

i = 1; Ntrial = 1;
while i < size(datamatrix,1)
   nexti = i + Nspikes(Ntrial);
   DATA(Ntrial+1,1:Nspikes(Ntrial) + Nparams+1) = [...		% first Ntrial reserved
         datamatrix(i,4:2+Nparams) ...		% params
         datamatrix(i,3) ...					% repnum
         Nspikes(Ntrial) ...						% Nspikes
         datamatrix(i:nexti-1,2)'];				% spiketimes
   i = nexti;
	Ntrial = Ntrial+1;   
end

switch nargin
    case 5
        eval(['save ' FN.data_path FN.data '.mat FN params DATA param3 param4 param5']);
    case 4
        eval(['save ' FN.data_path FN.data '.mat FN params DATA param3 param4']);
    case 3   
        eval(['save ' FN.data_path FN.data '.mat FN params DATA param3']);
    case 2
        eval(['save ' FN.data_path FN.data '.mat FN params DATA']);
end

