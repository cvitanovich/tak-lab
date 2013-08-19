%setInfo_itd

% increment test number
if get(H.itd_buildOnly,'value') | exist1('H.composite_fig')
    XStimParams.buildOnly = 1;
else
    XStimParams.buildOnly = 0;
end   
set(H.itd_buildOnly,'value',XStimParams.buildOnly);

if(exist1('H.itdfig') & get(H.itd_inc_testnum,'Value'))
	XStimParams.testnum = XStimParams.testnum +1;
	set(H.testnum, 'String', num2str(XStimParams.testnum))
   update_dataFN;
   set(H.itd_inc_testnum,'Value',0);
end

%setITDstep
if ~get(H.aritrary_itds,'value')
    low_itd  = str2num(get(H.lowitd,'String'));
    high_itd = str2num(get(H.highitd,'String'));
    if(high_itd > 250)
        high_itd = 250;
        set(H.highitd,'String',num2str(high_itd));
    end
    if(low_itd < -250)
        low_itd = -250;
        set(H.lowitd,'String',num2str(low_itd));
    end
    XStimParams.hiitd = high_itd;
    XStimParams.loitd = low_itd;
    XStimParams.numitds = str2num(get(H.numitds,'String'));   
    set(H.stepitd,'String',num2str(round((high_itd-low_itd)/(XStimParams.numitds-1))));
    XStimParams.itds = round(XStimParams.loitd:...
        (XStimParams.hiitd - XStimParams.loitd)/(XStimParams.numitds-1):...
        XStimParams.hiitd);
    set(H.itdvector,'string',num2str(XStimParams.itds));
else
    XStimParams.itds = str2num(get(H.itdvector,'string'));    
    XStimParams.numitds = length(XStimParams.itds);
    set(H.lowitd,'String','');
    set(H.highitd,'String','');
    set(H.stepitd,'String','');
end

set(H.numitds,'string',num2str(XStimParams.numitds));

% check out modulation parameters for stim
XStimParams.mod_depth(1) = str2num(get(H.itdtestmod_depth,'String'));
mod_type = get(H.itdtestmod_type,'String');
mod_num = get(H.itdtestmod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Tone'
        set(H.itdtestmod_txt,'String', 'Freq:');
        set(H.itdtestmod_freq,'Visible','on');
        set(H.itdtestmod_phase,'Visible','on');
        set(H.itdtestmod_txtA,'Visible','on');
        set(H.itdtestmod_depth,'Visible','on');
        set(H.itdtestmod_txtB,'Visible','on');
        set(H.itdtestmod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.itdtestmod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.itdtestmod_phase,'String'));
    case 'File'
        set(H.itdtestmod_pb,'Visible','on');
        set(H.itdtestmod_freq,'Visible','off');
        set(H.itdtestmod_txtA,'Visible','on');
        set(H.itdtestmod_depth,'Visible','on');
        set(H.itdtestmod_phase,'Visible','off');
        set(H.itdtestmod_txtB,'Visible','off');
        if get(H.itdtestmod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.itdtestmod_pb,'Value',0);
        end
        set(H.itdtestmod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.itdtestmod_txt,'String', 'CutOff Freq:');
        set(H.itdtestmod_freq,'Visible','on');
        set(H.itdtestmod_txtA,'Visible','on');
        set(H.itdtestmod_depth,'Visible','on');
        set(H.itdtestmod_phase,'Visible','off');
        set(H.itdtestmod_pb,'Visible','off');
        set(H.itdtestmod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.itdtestmod_freq,'String'));
        
    case 'None'
        set(H.itdtestmod_txt,'String', 'no mod  ');
        set(H.itdtestmod_freq,'Visible','off');
        set(H.itdtestmod_phase,'Visible','off');
        set(H.itdtestmod_pb,'Visible','off');
        set(H.itdtestmod_txtB,'Visible','off');
        set(H.itdtestmod_txtA,'Visible','off');
        set(H.itdtestmod_depth,'Visible','off');
        
    otherwise
end



% save current parameter values
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim