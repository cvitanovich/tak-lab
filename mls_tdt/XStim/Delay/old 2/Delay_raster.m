function [] = Delay_raster(params,datamatrix,seq, reps)
% Delay_raster 

global FN
global H

figure(H.delay_raster);

repnum=datamatrix(1,3);
%disp(repnum);

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

spikes = [];
rows = [];

i = 1; Ntrial = 1;
while i < size(datamatrix,1)
   nexti = i + Nspikes(Ntrial);
   DATA(Ntrial+1,1:Nspikes(Ntrial) + Nparams+1) = [...		% first Ntrial reserved
        datamatrix(i,4:2+Nparams) ...		% params
        datamatrix(i,3) ...				    % repnum
        Nspikes(Ntrial) ...				    % Nspikes
        datamatrix(i:nexti-1,2)'];			% spiketimes

        spikes=[spikes, [datamatrix(i:nexti-1,2)']];
        rows=[rows, [zeros(1,Nspikes(Ntrial)) + seq(Ntrial)] + ((0.5*(repnum/reps))-0.25)];
   i = nexti;
   Ntrial = Ntrial+1;   
end

hold on;
% spikes2 = [];
% rows2 = [];
% plot(spikes2, rows2,...
%     'Marker', '.',...
%     'LineStyle', 'none',...
%     'Color', [0 0 0]);
plot(spikes, rows,...
    'Marker', '*',...
    'MarkerFaceColor', [0 0 0],...
    'MarkerSize', 2,...
    'LineStyle', 'none',...
    'Color', [1 1 1]);
hold off;