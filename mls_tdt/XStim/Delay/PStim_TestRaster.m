function PStim_TestRaster(TrialData, sortColumn, PStimG)
    % PStim_FullRaster
    
    if nargin < 2
        sortColumn = [];
    end
    
   % global PStimG
   
    binWidth = 1;
    numTrials = size(TrialData,1);
    
    TrialData = sortrows(TrialData, sortColumn);
    %assignin('base', 'TrialData', TrialData);
    
    Spikes = TrialData(:,22:size(TrialData,2));
    Spikes = Spikes';
    [i] = find(Spikes == 0);
    Spikes(i) = NaN;
    Spikes = Spikes-PStimG.PreStimDelay;
    %assignin('base', 'Spikes', Spikes);
    
    if ~isempty(findobj('Tag','Test Raster'));
        close('Test Raster');
        clear TestRaster
    end
    
    Bpos=200; % bottom
    Lpos=440; % left
    Wpos=400; % width
    Hpos=800; % height
    
    axesColor = [0.7 0.7 0.7];
    sortLineColor = [0.5 0.5 0.5];
    
    TestRaster.fig = figure('Units','pixels',...
        'Position',[Lpos Bpos Wpos Hpos],...
        'Name','Test Raster',...
        'Tag','Test Raster',...
        'NumberTitle','off',...
        'MenuBar','figure',...
        'Color',[0 0 0],...
        'ToolBar','figure');
       
    % secondary sort order text
    if length(sortColumn) >= 2
        sr = TrialData(:,sortColumn(2));
        sr = sort(sr);
        sr = sr(~isnan(sr));
        sr = sr';
        sr = [sr NaN];
        sr = sr';
        srD = diff(sr);
        srL = sr(find(srD));
        if length(srL) < 20
            sortColumntxt = uicontrol(...
                'Style','text',...
               'Units','pixels',...
               'Position',[Wpos-25 -(Hpos/4) 25 Hpos],...
               'BackgroundColor',[0 0 0],...
               'ForegroundColor', axesColor,...
               'FontWeight', 'normal',...
               'FontSize', 8,...
               'HorizontalAlignment', 'right',...
               'String', num2str(srL));
        end
    end      
    
    % Primary sort break lines
    if length(sortColumn) >= 1
        sr = TrialData(:,sortColumn(1));
        sr = sr(~isnan(sr));
        sr = sr';
        sr = [sr NaN];
        sr = sr';
        srD = diff(sr);
        srL = sr(find(srD));
        sr = find(srD);
        
        for b=1:length(sr)
            line([-PStimG.PreStimDelay PStimG.StimDur+PStimG.PostStimDelay/2],[sr(b)-0.5 sr(b)-0.5],...
                'Marker', 'none',...
                'Color', sortLineColor,...
                'LineStyle', '-');

            if length(srL) < 20
                text(PStimG.StimDur+PStimG.PostStimDelay/2+2, sr(b)-0.5, num2str(srL(b)),...
                    'FontSize', 8,...
                    'Color', axesColor);
            end
        end
    end

     
    % stim lines
    line([0 0],[0 size(Spikes,2)],...
        'Marker', 'none',...
        'Color', [0.7 0 0],...
        'LineStyle', '-');
    line([PStimG.StimDur PStimG.StimDur],[0 size(Spikes,2)],...
        'Marker', 'none',...
        'Color', [0.7 0 0],...
        'LineStyle', '-');
    
    % Plot spikes
    r = 1:size(Spikes,2);
    r = repmat(r,size(Spikes,1),1);
    TestRaster.plot = line(Spikes, r,...
        'Marker', '.',...
        'MarkerFaceColor', [1 1 1],...
        'MarkerEdgeColor', [1 1 1],...
        'MarkerSize', 2,...
        'LineStyle', 'none');
    
        
    set(gca,'Box','off',...
        'XLim', [-PStimG.PreStimDelay/2 PStimG.StimDur+PStimG.PostStimDelay/2],...
        'YLim', [0 size(Spikes,2)],...
        'XColor', axesColor,...
        'YColor', axesColor,...
        'Color', [0 0 0]);
    
    xlabel('time (ms)', 'Color', axesColor);
    ylabel('trial#', 'Color', axesColor);

    
    clear Bpos Lpos Wpos Hpos ans Spikes