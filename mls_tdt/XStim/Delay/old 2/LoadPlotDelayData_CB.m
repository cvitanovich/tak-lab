% LoadPlotDelayData_CB
% CallBack for LoadPlotDelayData

switch lower(gcbo)
    case findobj('Tag', 'DData_file')
        LoadPlotDelayData;
        PlotDData;
        
    case findobj('Tag', 'DData_Plot')
        PlotDData;
    case findobj('Tag', 'DData_NamePSTH')
        %disp(TrialStr);
        %disp([TrialStr ' = spikes_N;']);
        %eval([TrialStr ' = spikes_N;']);
        
    otherwise
        warning('no object found');
end