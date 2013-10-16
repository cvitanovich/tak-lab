function PStim_SRaster(PStimG, RasterRows,Row,Times)
    % PStim_PlotSpikes
    % plots spikes as they are collected
    
    global SRaster
  
    RasterExists = ~isempty(findobj('Tag','Raster'));
     
    if Row==1
        if RasterExists
            figure(SRaster.fig);
            cla;
            %close('Raster');
        else
            Bpos=570; % bottom
            Lpos=440; % left
            Wpos=400; % width
            Hpos=300; % height
            
            SRaster.fig = figure('Units','pixels',...
              'Position',[Lpos Bpos Wpos Hpos],...
              'Name','Raster',...
              'Tag','Raster',...
              'NumberTitle','off',...
              'MenuBar','none',...
              'Color',[0 0 0],...
              'ToolBar','none');
          
              RasterExists=1;
        end
    end
    
    % plot only when front window?
%      if SRaster.fig ~= gcf
%          return;
%      end
    
    if RasterExists
        if Row ~= 1
            figure(SRaster.fig);
            hold on;
        
            % stim lines
            line([0 0],[0 RasterRows+1],...
                'Marker', 'none',...
                'Color', [1 0 0],...
                'LineStyle', '-');
            line([PStimG.StimDur PStimG.StimDur],[0 RasterRows+1],...
                'Marker', 'none',...
                'Color', [1 0 0],...
                'LineStyle', '-');
        end

        % raster
        line(Times-PStimG.PreStimDelay,Row,...
            'Marker', '.',...
            'MarkerFaceColor', [0.6 0.6 1],...
            'MarkerEdgeColor', [0.6 0.6 1],...
            'MarkerSize', 5,...
            'LineStyle', 'none');
        
        if Row==1
            set(gca,'Box','off',...
                'YLim', [0 RasterRows+1],...
                'XLim', [-PStimG.PreStimDelay/2 PStimG.StimDur+PStimG.PostStimDelay/2],...
                'XColor', [0.7 0.7 0.7],...
                'YTickLabel', '',...
                'Color', [0 0 0]);
        end
        
        
    end
    
    clear Bpos Lpos Wpos Hpos