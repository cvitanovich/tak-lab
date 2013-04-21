function PlotDData()
% plot data across reps 

global DData;
global PSTH;
global Rast;

DATA = evalin('base','DATA');
reps = str2num(get(DData.reps,'String'));
params=(DATA(1,1));
if(params~=8)
    warning('function not compatible with data collected');
    return;
end

p1=str2num(get(DData.elev,'String'));
p2=str2num(get(DData.azim,'String'));

p3=str2num(get(DData.delay,'String'));
p4=str2num(get(DData.onset,'String'));
p5=str2num(get(DData.offset,'String'));
p6=str2num(get(DData.LagExt,'String'));
p7=str2num(get(DData.DelayMod,'String'));

spikes = [];
repnum = [];
numtrials = size(DATA,1);
for tr=2:numtrials  % first row is num params
    if DATA(tr,1)==p1 & ...
        DATA(tr,2)==p2 & ...
        DATA(tr,3)==p3 & ...
        DATA(tr,4)==p4 & ...
        DATA(tr,5)==p5 & ...
        DATA(tr,6)==p6 & ...
        DATA(tr,7)==p7

        spikes = [spikes DATA(tr, 10:10+DATA(tr,9)-1)];
        repnum = [repnum (DATA(tr, 10:10+DATA(tr,9)-1)*0)+ DATA(tr,8)];
    end
end

if(size(spikes,2))
    set(DData.FoundTrial, 'Visible', 'off');
    TrialStr = ['T_' num2str(p1) '_' num2str(p2) '_' num2str(p3) '_' ...
        num2str(p4) '_' num2str(p5) '_' num2str(p6) '_' num2str(p7)];
    assignin('base','TrialStr',TrialStr);
else
    set(DData.FoundTrial, 'Visible', 'on');
end


% spikes_all=[];
% repnum_all=[];
% sn=1;
% for s=1:size(spikes,2)
%    if spikes(s)~=0
%        spikes_all(sn)=spikes(s);
%        repnum_all(sn)=repnum(s);
%        sn=sn+1;
%    end
% end
spikes_all=spikes;
repnum_all=repnum;

max_spike_time = max(max(spikes));

% plot raster
if( get(DData.raster, 'Value') == 1 );
    if(findobj('Tag','Trial Raster'))
        close('Trial Raster');
    end
    Rast.fig = figure('Units','pixels',...
          'Position',[700 50 400 200],...
          'Tag', 'Trial Raster',...
          'Name','Trial Raster',...
          'NumberTitle','off',...
          'Color',[1 1 1]);
        Rast.axes = axes;
        figure(Rast.fig);
        ylabel('trial'); xlabel('time (ms)');
        ylim([0 numtrials+1]);
        set(Rast.axes, 'YLim', [0 reps+0.5]);

    hold on;
    Rast.spikes = plot(spikes_all, repnum_all,...
        'Marker', '*',...
        'MarkerFaceColor', [0 0 0],...
        'MarkerSize', 2,...
        'LineStyle', 'none',...
        'Color', [0 0 0]);
end


% plot PSTH
if( get(DData.PSTH, 'Value') == 1 );
    binWidth=5; %ms
    max_spike_time=ceil(max_spike_time/binWidth)*binWidth;
    
    spikes_bin = [0:binWidth:max_spike_time];
    spikes_N = zeros(1,size(spikes_bin,2));

    for b=1:size(spikes_bin,2)
        for s=1:size(spikes_all,2)
           if spikes_all(s) >= spikes_bin(b) & spikes_all(s) < spikes_bin(b+1)
               spikes_N(b)=spikes_N(b)+1;
           end
        end
    end
    spikes_bin=spikes_bin+(binWidth/2);

   assignin('base','spikes_bin',spikes_bin');
   assignin('base','spikes_N',spikes_N');

    if( get(DData.holdPSTH, 'Value') & findobj('Tag','PSTH'));
        figure(PSTH.fig);
        hold on;
    else
        if(findobj('Tag','PSTH'))
            close('PSTH');
        end
        PSTH.fig = figure('Units','pixels',...
          'Position',[700 350 400 200],...
          'Tag', 'PSTH',...
          'Name','PSTH',...
          'NumberTitle','off',...
          'Color',[1 1 1]);
        
    end
    if(0)
        Rast.spikes = bar(spikes_bin, spikes_N, 1);
    else
        Rast.spikes = plot(spikes_bin, spikes_N,...
            'Marker', 'none',...
            'MarkerFaceColor', [0 0 0],...
            'MarkerSize', 2,...
            'LineStyle', '-',...
            'Color', [0 0 0]);
        
        ylabel('# spikes'); xlabel('time (ms)');
    end
        set(GCA, 'Box', 'off');
        
    if( get(DData.holdPSTH, 'Value') );
       hold off;
    end

end