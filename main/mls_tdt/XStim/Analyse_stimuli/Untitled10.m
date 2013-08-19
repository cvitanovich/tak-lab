        figure
        m = min1([20*log10(mean(Vstr_ABL_20,2)) 20*log10(mean(Vstr_ABL_55,2)) 20*log10(mean(Vstr_ABL_75,2))])
        M = max1([20*log10(mean(Vstr_ABL_20,2)) 20*log10(mean(Vstr_ABL_55,2)) 20*log10(mean(Vstr_ABL_75,2))])
        % 20 Hz
        h = axes('position',[.24 (.01 + (4)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(20*log10(mean(Vstr_ABL_20,2)), dir, 1);
        caxis([m M])
        % 55 Hz
        h = axes('position',[.44 (.01 + (4)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(20*log10(mean(Vstr_ABL_55,2)), dir, 1);
        caxis([m M])
        % 75 Hz
        h = axes('position',[.64 (.01 + (4)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(20*log10(mean(Vstr_ABL_75,2)), dir, 1);
        caxis([m M])
       
        m = min1([20*log10(mean(Vstr_dILD_20,2)) 20*log10(mean(Vstr_dILD_55,2)) 20*log10(mean(Vstr_dILD_75,2))])
        M = max1([20*log10(mean(Vstr_dILD_20,2)) 20*log10(mean(Vstr_dILD_55,2)) 20*log10(mean(Vstr_dILD_75,2))])
        % 20 Hz
        h = axes('position',[.24 (.01 + (3)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(20*log10(mean(Vstr_dILD_20,2)), dir, 1);
        caxis([m M])
        % 55 Hz
        h = axes('position',[.44 (.01 + (3)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(20*log10(mean(Vstr_dILD_55,2)), dir, 1);
        caxis([m M])
        % 75 Hz
        h = axes('position',[.64 (.01 + (3)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(20*log10(mean(Vstr_dILD_75,2)), dir, 1);
        caxis([m M])

        m = min1([20*log10(mean(Vstr_dIPD_20,2)) 20*log10(mean(Vstr_dIPD_55,2)) 20*log10(mean(Vstr_dIPD_75,2))])
        M = max1([20*log10(mean(Vstr_dIPD_20,2)) 20*log10(mean(Vstr_dIPD_55,2)) 20*log10(mean(Vstr_dIPD_75,2))])
        % 20 Hz
        h = axes('position',[.24 (.01 + (2)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(20*log10(mean(Vstr_dIPD_20,2)), dir, 1);
        caxis([m M])
        % 55 Hz
        h = axes('position',[.44 (.01 + (2)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(20*log10(mean(Vstr_dIPD_55,2)), dir, 1);
        caxis([m M])
        % 75 Hz
        h = axes('position',[.64 (.01 + (2)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(20*log10(mean(Vstr_dIPD_75,2)), dir, 1);
        caxis([m M])

        
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        figure
        % 20 Hz
        h = axes('position',[.24 (.01 + (2)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        temp = (20*log10(mean(Vstr_dIPD_20,2))) + (20*log10(mean(Vstr_dILD_20,2)));
        plot_diam_2axis(temp, dir, 1);
        % 55 Hz
        h = axes('position',[.44 (.01 + (2)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        temp = (20*log10(mean(Vstr_dIPD_55,2))) + (20*log10(mean(Vstr_dILD_55,2)));
        plot_diam_2axis(temp, dir, 1);
        % 75 Hz
        h = axes('position',[.64 (.01 + (2)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        temp = (20*log10(mean(Vstr_dIPD_75,2))) + (20*log10(mean(Vstr_dILD_75,2)));
        plot_diam_2axis(temp, dir, 1);
%%%%%%%%%%%%%%%%%%%%%%%%% Zscores

        figure
        temp = [mean(Z_ABL_20,2) mean(Z_ABL_55,2) mean(Z_ABL_75,2) ...
            mean(Z_dILD_20,2) mean(Z_dILD_55,2) mean(Z_dILD_75,2) ...
            mean(Z_dIPD_20,2) mean(Z_dIPD_55,2) mean(Z_dIPD_75,2)];
        m = min1(temp)
        M = max1(temp)
        % 20 Hz
        h = axes('position',[.24 (.01 + (4)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(mean(Z_ABL_20,2), dir, 1);
        caxis([m M])
        % 55 Hz
        h = axes('position',[.44 (.01 + (4)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(mean(Z_ABL_55,2), dir, 1);
        caxis([m M])
        % 75 Hz
        h = axes('position',[.64 (.01 + (4)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(mean(Z_ABL_75,2), dir, 1);
        caxis([m M])
       
        % 20 Hz
        h = axes('position',[.24 (.01 + (3)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(mean(Z_dILD_20,2), dir, 1);
        caxis([m M])
        % 55 Hz
        h = axes('position',[.44 (.01 + (3)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(mean(Z_dILD_55,2), dir, 1);
        caxis([m M])
        % 75 Hz
        h = axes('position',[.64 (.01 + (3)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(mean(Z_dILD_75,2), dir, 1);
        caxis([m M])

        % 20 Hz
        h = axes('position',[.24 (.01 + (2)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(mean(Z_dIPD_20,2), dir, 1);
        caxis([m M])
        % 55 Hz
        h = axes('position',[.44 (.01 + (2)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(mean(Z_dIPD_55,2), dir, 1);
        caxis([m M])
        % 75 Hz
        h = axes('position',[.64 (.01 + (2)*.20) .21 .21]);
        set(h,'XTickLabel',[])
        set(h,'YTickLabel',[])
        plot_diam_2axis(mean(Z_dIPD_75,2), dir, 1);
        caxis([m M])

