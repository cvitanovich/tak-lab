function Delay_PlotData(TrialData, PStimG)
    % Delay_PlotData
    % calculates interactions indicies across trials
    
    %global PStimG

    %Delays = PStimG.StimDelays;
    Delays = TrialData(:,9);
    Delays = sort(Delays);
    d = diff(Delays);
    Delays = Delays(find(d));
    Delays = Delays(~isnan(Delays));
        
    Mods = TrialData(:,13);
    Mods = sort(Mods);
    d = diff(Mods);
    Mods = Mods(find(d));
    Mods = Mods(~isnan(Mods));

    data = NaN*zeros(length(Delays),length(Mods));
    dataCnt = zeros(length(Delays),length(Mods))*0+0.0000000001;
        
    NaNSum = 0;
    NaNCnt = 0;
    
    FirstSpikeTime = PStimG.PreStimDelay;
    LastSpikeTime = PStimG.PreStimDelay + PStimG.StimDur + 50; % adding 50 to include offset response...
    
    for tr = 1:size(TrialData,1)
        if 1 % windowed spikes
            SpikesInWindow = find(TrialData(tr,22:size(TrialData,2)) > FirstSpikeTime & TrialData(tr,22:size(TrialData,2)) < LastSpikeTime );
            SpikesInWindow = length(SpikesInWindow);
            
            % in case data were converted from XStim and only spontaneous spike count was recorded?
            if isnan(TrialData(tr,9)) & SpikesInWindow == 0 & TrialData(tr,21) > 0
                temp = find(TrialData(tr,22:size(TrialData,2)) > 0 );
                temp = length(temp);
                if temp == 0 % still zero
                    SpikesInWindow = TrialData(tr,21);
                end
            end
            
        else % all spikes 
            SpikesInWindow = TrialData(tr,21);
        end
        
        Di = find(Delays == TrialData(tr,9));
        Mi = find(Mods == TrialData(tr,13));
        
        if ~isnan(TrialData(tr,9)) & ~isnan(TrialData(tr,13)) & ~isnan(data(Di,Mi)) % not looking at spontaneous rate
            data(Di,Mi) = data(Di,Mi) + SpikesInWindow;
            dataCnt(Di,Mi) = dataCnt(Di,Mi) + 1;
        elseif ~isnan(TrialData(tr,9)) & ~isnan(TrialData(tr,13)) & isnan(data(Di,Mi))
            data(Di,Mi) =  SpikesInWindow;
            dataCnt(Di,Mi) =  1;
        else
            if ~isnan(TrialData(tr,21))
                NaNSum = NaNSum + SpikesInWindow;
                NaNCnt=NaNCnt+1;
            end
        end
        
    end
    
    data = data./dataCnt;
        
    % subtracts spontaneous rate
    baseline = NaNSum/NaNCnt;
    data = data - baseline;
     
    Bpos=50; % bottom
    Lpos=850; % left
    Wpos=400; % width
    Hpos=300; % height
    
    if ~isempty(findobj('Tag','Delay plot'));
        close('Delay plot');
    end

    figColor = [0 0 0];
    axesColor = [0.8 0.8 0.8];
    
    DelayInt.fig = figure('Units','pixels',...
      'Position',[Lpos Bpos Wpos Hpos],...
      'Name','Delay plot',...
      'Tag','Delay plot',...
      'NumberTitle','off',...
      'MenuBar','figure',...
      'Color',figColor,...
      'ToolBar','auto');
    
   
    cmap = jet;
    set(0,'DefaultAxesColorOrder',[cmap(5,:);cmap(21,:);cmap(43,:);cmap(59,:)],...
      'DefaultAxesLineStyleOrder','-|:|--')
    hold on;
     
    DelayInt.mods = plot(Delays, data,...
        'Marker', '.',...
        'MarkerSize', 5,...
        'LineWidth', 0.5);
    
    set(gca,'Box','off',...
        'XColor', axesColor,...
        'YColor', axesColor,...
        'Color', figColor);
            
    LabelStr = num2str(Mods(1));
    if length(Mods)>1
        for n = 2:length(Mods)
            LabelStr= strvcat(LabelStr,num2str(Mods(n)));
        end
    end
   
   % DelayInt.legend = legend(LabelStr, 'Location', 'NorthOutside','Orientation','horizontal');
    DelayInt.legend = legend(LabelStr);  % above calls not supported in Matlab < 7.0???
    legend('boxoff'), legend(gca,'boxoff');
    set(findobj('type', 'text'), 'color', axesColor);
    
    hold off;
        
    clear ans Bpos Lpos Wpos Hpos
    
    
    