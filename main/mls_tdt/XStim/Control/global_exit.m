% Global_exit
% script to exit from everything but Xstim control

% initialize secondary process
try(S232('getS2primary'));
catch
    if(S232('S2init', 0, 'INIT_PRIMARY', 20000) == 0)
        disp('Cannot initialize a primary process')
        return;
    end
end

% get locks
    if(S232('APlock', 100, 0) == 0)
        disp('Cannot acquire lock on AP2 Card')
        s232('S2close');
        return;
    end
    if(S232('XBlock', 100, 0) == 0)
        disp('Cannot acquire lock on X Bus')
        s232('APunlock', 0);
        s232('S2close');
        return;
    end

% close process, remove locks ...
if(S232('getAPlockstatus') & S232('getXBlockstatus'))
    S232('PA4mute',1);
    S232('PA4mute',2);
    S232('PD1stop',1);
    S232('PD1clear',1);
    S232('trash');
    S232('dropall');
    S232('S2close')              % close application and release AP2 and XBUS locks
end
m110dx( C_.STOP);

% reset flags on XStim control
set(H.ephonefile,'Enable','on');
set(H.ephoneuseit,'Enable','on');
set(H.locfile,'Enable','on');
set(H.locAZ,'Enable','on');
set(H.locEL,'Enable','on');
set(H.locuseit,'Enable','on');

% reset XStimParams for given test
GUI.UseLastLocations_flag =0;
XStimParams.locations = [];

switch XStimParams.test_type
    case 'Search for Cells'
        
    case 'ABL'
        XStimParams.loabl = str2num(get(H.lowabl,'String'));
        XStimParams.hiabl = str2num(get(H.highabl,'String'));
        XStimParams.numabls = str2num(get(H.numabls,'String'));
        XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
        XStimParams.curr_ILD = str2num(get(H.ILD,'String'));
        XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
    case 'ITD'
        XStimParams.loitd = str2num(get(H.lowitd,'String'));
        XStimParams.hiitd = str2num(get(H.highitd,'String'));
        XStimParams.numitds = str2num(get(H.numitds,'String'));
        XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
        XStimParams.curr_ILD = str2num(get(H.ILD,'String'));
        XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
    case 'FREQ'
        XStimParams.lofreq = str2num(get(H.lowfreq,'String'));
        XStimParams.hifreq = str2num(get(H.highfreq,'String'));
        XStimParams.numfreqs = str2num(get(H.numfreqs,'String'));
        XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
        XStimParams.curr_ILD = str2num(get(H.ILD,'String'));
        XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
        XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
    case 'ILD_FREQ'
        XStimParams.lofreq = str2num(get(H.lowfreq,'String'));
        XStimParams.hifreq = str2num(get(H.highfreq,'String'));
        XStimParams.numfreqs = str2num(get(H.numfreqs,'String'));
        XStimParams.loild = str2num(get(H.lowild,'String'));
        XStimParams.hiild = str2num(get(H.highild,'String'));
        XStimParams.numilds = str2num(get(H.numilds,'String'));
        XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
        XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
        XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
    case 'ABI_FREQ'
        
    case 'Space'
        XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
        XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
        XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
    case 'Space2'
        XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
        XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
        XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
    case 'Space3'
        XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
        XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
    case '2Source'
        XStimParams.bandwidth(1) = str2num(get(H.bandwidth1,'String'));
        XStimParams.bandlimits(1) = str2num(get(H.bandlimits1,'String'));
        XStimParams.bandlimits(2) = str2num(get(H.bandlimits2,'String'));
        XStimParams.bandwidth(2) = str2num(get(H.bandwidth2,'String'));
        XStimParams.bandlimits(3) = str2num(get(H.bandlimits3,'String'));
        XStimParams.bandlimits(4) = str2num(get(H.bandlimits4,'String'));
        XStimParams.bandwidth(3) = str2num(get(H.bandwidth3,'String'));
        XStimParams.bandlimits(5) = str2num(get(H.bandlimits5,'String'));
        XStimParams.bandlimits(6) = str2num(get(H.bandlimits6,'String'));
        XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
        XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
        XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
    case 'AltIR'
        XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
        XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
        XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
    case 'Delay'
        XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
    case 'Multisource'
        
    case 'AM'
        XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
        XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
        XStimParams.test_ISI = str2num(get(H.ISI,'String'));
        XStimParams.numreps = str2num(get(H.numreps,'String'));
        
end

% close figures
if(exist1('H.rasterfig'))   close(H.rasterfig); H.rasterfig =[];    end
if(exist1('H.finalspikeratefig'))   close(H.finalspikeratefig); H.finalspikeratefig =[]; end
if exist1('H.pickerfig')    close(H.pickerfig); H.pickerfig = [];	end
if(exist1('H.ablfig'))      close(H.ablfig); H.ablfig =[];          end
if(exist1('H.searchfig'))   close (H.searchfig); H.searchfig =[]; H.psdaxes=[]; end
if(exist1('H.spikefig'))    close(H.spikefig); H.spikefig=[]; H.spikeaxes =[]; end
if(exist1('H.itdfig'))      close(H.itdfig); H.itdfig =[];          end
if(exist1('H.freqtestfig')) close(H.freqtestfig); H.freqtestfig =[]; end
if(exist1('H.ildfreqfig'))  close(H.ildfreqfig); H.ildfreqfig =[];  end
if(exist1('H.spacefig'))    close(H.spacefig); H.spacefig =[];      end
if exist1('H.wts_fig')      close(H.wts_fig);  H.wts_fig = [];      end
if exist1('H.IAresponse_fig')   close(H.IAresponse_fig); H.IAresponse_fig=[];   end
if(exist1('H.spacefig3'))   close(H.spacefig3); H.spacefig3 =[];    end
if(exist1('H.Two_sourcefig')) close(H.Two_sourcefig); H.Two_sourcefig =[]; end
if(exist1('H.altIRfig'))    close(H.altIRfig); H.altIRfig =[];      end
if(exist1('H.Delayfig'))    close(H.Delayfig); H.Delayfig =[];      end
if(exist1('H.multi_sourcefig')) close(H.multi_sourcefig); H.multi_sourcefig =[]; end
if(exist1('H.setparamsfig')) close(H.setparamsfig); H.setparamsfig =[]; end
if(exist1('H.AM_fig'))      close(H.AM_fig); H.AM_fig =[];          end

XStimParams.test_type = 'None';
set(H.test_type,'Value',1);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ' FN.current_path 'FN_current FN;'])

clear s* t*
clear GLOBAL H
eval(['load ' FN.current_path 'H_xstim;'])