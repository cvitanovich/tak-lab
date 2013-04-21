% lookatData

% Make the stimuli and then plot their spike probabilities

%global DataPlot

if 1 % make the stimuli (or else just make the figure)
    % allocate space
    leadProb = zeros(1,size(data,1));
    lagProb = zeros(1,size(data,1));
    
    for tr = 1:size(data,1)
        [Noise1, Noise2, Env1, Env2] = filteredstimuliAM(data, tr, 0);
        [leadP, lagP] = SpikeProbabilitiesAM(Env1, Env2, 0, Noise1, Noise2);
        leadProb(tr) = leadP;
        lagProb(tr) = lagP; 
    end
    %beep
end

plotType = 1; % 1=lead-lag, 2=lead, 3=lag
for plotType = 1:3
    
    % plot the spike probability distributions
    figColor = [1 1 1];
    if isempty(findobj('Tag',['Data plot ' num2str(plotType)]))
    %    DataPlot.fig = figure('Units','pixels',...
        figure('Units','pixels',...
          'Position',[400 500 72*6 72*3],...
          'Name',['Data plot ' num2str(plotType)],...
          'MenuBar','figure',...
          'ToolBar','figure',...
          'Tag',['Data plot ' num2str(plotType)],...
          'NumberTitle','off',...
          'Color',figColor);
    else
        figure(findobj('Tag',['Data plot ' num2str(plotType)])); %DataPlot.fig);
        clf;
    end

    cmap = colormap(hsv(5));
    hold on
    for p = 0:length(freqs)-1
        first = p*length(delays)+1;
        last = first-1+length(delays);
        if plotType == 1
            line(delays, leadProb(first:last)-lagProb(first:last), 'color', cmap(p+1,:),'Marker','.'); % lead - lag
        elseif plotType == 2
            line(delays, leadProb(first:last), 'color', cmap(p+1,:),'Marker','.'); % lead
        elseif plotType == 3
            line(delays, lagProb(first:last), 'color', cmap(p+1,:),'Marker','.'); % lag
        end
    end
    DataPlot.legend = legend('0.25 kHz','0.5','1','2','6'); % check for appropriate values...
    legend(DataPlot.legend,'boxoff');
    set(gca,'XTick', 0:1:max(delays));
    set(gca,'XLim', [-0.5 max(delays)]);
    switch plotType
        case 1
            ylabel('lead - lag response')
        case 2
            ylabel('lead response')
        case 3
            ylabel('lag response');
    end
end
clear figColor plotType
clear cmap tr p first last
