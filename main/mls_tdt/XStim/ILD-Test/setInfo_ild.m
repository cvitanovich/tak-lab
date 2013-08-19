%setInfo_ild

% increment test number
if(exist1('H.ildfig') & get(H.ild_inc_testnum,'Value'))
	XStimParams.testnum = XStimParams.testnum +1;
	set(H.testnum, 'String', num2str(XStimParams.testnum))
   update_dataFN;
   set(H.inc_testnum,'Value',0);
end

if get(H.ild_buildOnly,'value') | exist1('H.composite_fig')
    XStimParams.buildOnly = 1;
else
    XStimParams.buildOnly = 0;
end   
set(H.ild_buildOnly,'value',XStimParams.buildOnly);

%SetILDStep
if ~get(H.aritrary_ilds,'value')
    XStimParams.loild  = str2num(get(H.lowild,'String'));
    XStimParams.hiild = str2num(get(H.highild,'String'));
    XStimParams.numilds = str2num(get(H.numilds,'String'));
    set(H.stepild,'String',num2str(round((XStimParams.hiild-XStimParams.loild)/(XStimParams.numilds-1))));
    XStimParams.ilds = round(XStimParams.loild:...
        (XStimParams.hiild - XStimParams.loild)/(XStimParams.numilds-1):...
        XStimParams.hiild);
    set(H.ildvector,'string',num2str(XStimParams.ilds));
else
    XStimParams.ilds = str2num(get(H.ildvector,'string'));    
    XStimParams.numilds = length(XStimParams.ilds);
    set(H.lowild,'String','');
    set(H.highild,'String','');
    set(H.stepild,'String','');
end

set(H.numilds,'string',num2str(XStimParams.numilds));

XStimParams.useILDshape = get(H.useildshape,'value');
XStimParams.numreps = str2num(get(H.ild_numreps,'string'));

% check out modulation parameters for stim
XStimParams.mod_depth(1) = str2num(get(H.ildtestmod_depth,'String'));
mod_type = get(H.ildtestmod_type,'String');
mod_num = get(H.ildtestmod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Tone'
        set(H.ildtestmod_txt,'String', 'Freq:');
        set(H.ildtestmod_freq,'Visible','on');
        set(H.ildtestmod_phase,'Visible','on');
        set(H.ildtestmod_txtA,'Visible','on');
        set(H.ildtestmod_depth,'Visible','on');
        set(H.ildtestmod_txtB,'Visible','on');
        set(H.ildtestmod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.ildtestmod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.ildtestmod_phase,'String'));
    case 'File'
        set(H.ildtestmod_pb,'Visible','on');
        set(H.ildtestmod_freq,'Visible','off');
        set(H.ildtestmod_txtA,'Visible','on');
        set(H.ildtestmod_depth,'Visible','on');
        set(H.ildtestmod_phase,'Visible','off');
        set(H.ildtestmod_txtB,'Visible','off');
        if get(H.ildtestmod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.ildtestmod_pb,'Value',0);
        end
        set(H.ildtestmod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.ildtestmod_txt,'String', 'CutOff Freq:');
        set(H.ildtestmod_freq,'Visible','on');
        set(H.ildtestmod_txtA,'Visible','on');
        set(H.ildtestmod_depth,'Visible','on');
        set(H.ildtestmod_phase,'Visible','off');
        set(H.ildtestmod_pb,'Visible','off');
        set(H.ildtestmod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.ildtestmod_freq,'String'));
        
    case 'None'
        set(H.ildtestmod_txt,'String', 'no mod  ');
        set(H.ildtestmod_freq,'Visible','off');
        set(H.ildtestmod_phase,'Visible','off');
        set(H.ildtestmod_pb,'Visible','off');
        set(H.ildtestmod_txtB,'Visible','off');
        set(H.ildtestmod_txtA,'Visible','off');
        set(H.ildtestmod_depth,'Visible','off');
        
    otherwise
end



% save current parameter values
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim