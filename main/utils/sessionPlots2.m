function sessionPlots(options)
% plots information about a PDR session in real time 
% session=structure with data for the session, 
% options=tells sessionPlots what action to take
% required to initialize:
% structure with these parameters:
% ntrials, min_yes, max_yes
% bufpts, decpts, isi, trials_to_show
% zoomval

% other required parameters:
% sound_onset, npts_totalplay, srate
% last_buffer
% stim_fs
global session
switch options
    
    case 'Initialize'
        
        % INIT PLOTS
        % handles

        session.hTxt=zeros(1,10); % text handles
        session.hSub=zeros(1,10); % subplot handles
        session.hMark=zeros(1,100);
        session.hMark2=zeros(1,100);
        session.hBox=0; % box marker handle
        session.ymax=32768; session.ymin=-32768; % for plotting trace
        
        figure(session.hFig); subplot(session.hTracePlot);
        set(gca,'XtickLabel','');
        title('PDR Trace');
        session.decpts=ceil(session.bufpts/2^session.dec_fact);
        session.trace_pts=session.decpts*session.isi*session.trials_to_show;
        session.trace_xes=1:session.trace_pts;
        % circular buffered trace
        session.trace_yes=zeros(1,session.trace_pts);
        session.trace_xes(end+1)=session.trace_pts; 
        session.trace_yes(end+1)=session.ymax;
        session.trace_xes(end+1)=0;
        session.trace_yes(end+1)=session.ymax;
        session.dec_fs = ceil(session.stim_fs/2^session.dec_fact); % decimated sample rate
        %session.dec_xes=session.trace_xes./session.dec_fs; % convert to time values
        session.hLine=line([session.trace_xes session.trace_xes],...
            [session.trace_yes session.trace_yes],'LineWidth',1,'Color','w');
        xlim([0 session.trace_pts]); ylim([session.ymin session.ymax]);
        % initialize variable for latest buffer
        session.last_buffer=zeros(1,session.decpts);
        session.test_flag=0;
        session.flag_list=[];
        session.test_flag_list=[];
        drawnow;

    case 'Update Trace Plot'
        figure(session.hFig); subplot(session.hTracePlot); hold on;
        % circular buffer updatate
        session.trace_yes(1:(session.trace_pts-session.decpts))=session.trace_yes((session.decpts+1):end-2);
        session.trace_yes((end-session.decpts-1):(end-2))=session.last_buffer;
        
        % update xes
        %session.dec_xes(1:end-2)=session.dec_xes(1:end-2)+session.decpts;
        %session.dec_xes(end-1)=session.dec_xes(end-2);
        %session.dec_xes=session.trace_xes./session.dec_fs; % convert to time values
        session.flag_list=session.flag_list-session.decpts;
        tmp=find(session.flag_list>=0); session.flag_list=session.flag_list(tmp);
        session.test_flag_list=session.test_flag_list-session.decpts;
        tmp=find(session.test_flag_list>=0); session.test_flag_list=session.test_flag_list(tmp);
        if session.test_flag==1
            session.flag_list=[session.trace_xes(end-2-session.decpts+ceil(session.sound_onset/2^session.dec_fact)) session.flag_list]; %[session.dec_xes(end-2-session.decpts+1+ceil(session.sound_onset/2^session.dec_fact)) session.flag_list];
        elseif session.test_flag==Inf
            session.test_flag_list=[session.trace_xes(end-2-session.decpts+ceil(session.sound_onset/2^session.dec_fact)) session.test_flag_list];%[session.dec_xes(end-2-session.decpts+1+ceil(session.sound_onset/2^session.dec_fact)) session.test_flag_list];
        end
        % update trace plot
        if session.ymax>session.zoomval*max(session.trace_yes(1:end-2)) || session.ymin<session.zoomval*min(session.trace_yes(1:end-2))
            session.ymax=max(session.trace_yes(1:end-2))+session.zoomval*abs(max(session.trace_yes(1:end-2))); session.ymin=min(session.trace_yes(1:end-2))-session.zoomval*abs(min(session.trace_yes(1:end-2)));
            ylim([session.ymin session.ymax]);
            session.trace_yes(end-1:end)=session.ymax*ones(1,2);
        end
            
        set(session.hLine,'xdata',[session.trace_xes session.trace_xes],'ydata',...
            [session.trace_yes session.trace_yes]);
        % plot sound onset indicators
        %start=session.dec_xes(1);
        %stop=session.dec_xes(end-2);
        start=session.trace_xes(1);
        stop=session.trace_xes(end-2);
        
        % marker lines
        chk=1; hTemp=0;
        while chk
            
            hTemp=findobj(gca,'Tag','marker');
            if hTemp>0
                delete(hTemp(1)); hTemp=0;
            else
                break;
            end
        end
        
        temp=find(session.flag_list>=start);
        for i=1:length(temp)
            session.hMark(i)=line([session.flag_list(temp(i)) session.flag_list(temp(i))],...
                [session.ymin session.ymax],'LineWidth',1,'Color','g','LineStyle','-','Tag','marker');
        end
        temp2=find(session.test_flag_list>=start);
        for i=1:length(temp2)
            session.hMark2(i)=line([session.test_flag_list(temp2(i)) session.test_flag_list(temp2(i))],...
                [session.ymin session.ymax],'LineWidth',1,'Color','r','LineStyle','-','Tag','marker');
        end

        drawnow; % only redraw plot when stimuli are played

end

return;