function [DATA, params] = read_data3_minus_spont_RASTER1(dataFN,dataPATH, repflag, locs)

% function [DATA, params] = read_data3_minus_spont_RASTER1(dataFN,dataPATH,repflag,locs);
%
% to read data created by record_data3 in Xstim testing module
%% data was saved as:
%	param1	param2	param3 ...	repnum	Nspikes	spiketimes...
% first line is a dummy with the number of params (including repnum) written to param1
% argin repflag: allows testing of one or more reps (e.g. [1 3]; default is all reps)
% argin: locs (2 x nLocs) el, az for locations to be plotted as raster
%%%
% modified version to plot raster of specific loci
% eliminated much of the original code
% must be a spatial test

if nargin < 2 | isempty(dataPATH)   dataPATH = 'e:\kip\datastor\';     end

eval(['load ' dataPATH dataFN]);
Nparams = DATA(1,1);

if nargin < 3 | isempty(repflag)    repflag = 1:params.numreps;    end       % use all reps by default
Nreps = length(repflag);
maxspikes = max1(DATA(:,Nparams+1));

if ~exist1('params.silence_lead')
    params.silence_lead = 100;
end



                    % adjust spiketimes for sound's convolution with HRTFs
                    %spiketimes = spiketimes - 4.1;
                    % adjust spiketimes for difference in MII and TDT clocks
                    %spiketimes = spiketimes .* 1.0016;

switch params.test_type
    case {'Space', 'Space3', '2Source', 'AltIR'}
        if nargin < 4  | isempty(locs)
            locs = params.locations;
        end
        ele = locs(1,:); 
        azi = locs(2,:);
        figure; hold on;
        goodtrial = 1;
        for trial = 2:size(DATA,1)
            ind = find(ele == DATA(trial,1) & azi == DATA(trial,2));
                rep = DATA(trial,Nparams);
                if any(ismember(rep,repflag))
                    ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > params.silence_lead );
                    ind2 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + params.curr_stimdur));
                    %plot((DATA(trial,Nparams+2:Nparams + 1 + maxspikes)-4.1)*1.0016,trial* ones(maxspikes,1),'.')
                    if ~isempty(ind)
                        plot((DATA(goodtrial,Nparams+1+setdiff(ind1,ind2))-4.1)*1.0016,trial* ones(length(ind1)-length(ind2),1),'b.')
                        goodtrial = goodtrial+1;
                    else
                        %plot((DATA(trial,Nparams+1+setdiff(ind1,ind2))-4.1)*1.0016,trial* ones(length(ind1)-length(ind2),1),'g.')
                    end          % end of repflag for inclusion of only chosen reps
            end             % end of ~isempty(ind)
        end         % end of trial
        
    otherwise
        disp('ERROR: read_data3_minus_spont_RASTER1 requires spatial data');
        return;
end
