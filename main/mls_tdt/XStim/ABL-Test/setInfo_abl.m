%setInfo_abl

if get(H.abl_buildOnly,'value') | exist1('H.composite_fig')
    XStimParams.buildOnly = 1;
else
    XStimParams.buildOnly = 0;
end   
set(H.abl_buildOnly,'value',XStimParams.buildOnly);

% increment test number
if(exist1('H.ablfig') & get(H.inc_testnum,'Value'))
	XStimParams.testnum = XStimParams.testnum +1;
	set(H.testnum, 'String', num2str(XStimParams.testnum))
   update_dataFN;
   set(H.inc_testnum,'Value',0);
end

%set ABL-step
if ~get(H.aritrary_abls,'value')
    low_abl  = str2num(get(H.lowabl,'String'));
    high_abl = str2num(get(H.highabl,'String'));
    if(high_abl > -20)
        high_abl = -20;
        set(H.highabl,'String',num2str(high_abl));
    end
    XStimParams.numabls = str2num(get(H.numabls,'String'));
    XStimParams.hiabl = high_abl;
    XStimParams.loabl = low_abl;
    XStimParams.numabls = str2num(get(H.numabls,'String'));   
    set(H.stepabl,'String',num2str(round((high_abl-low_abl)/(XStimParams.numabls-1))));
    XStimParams.abls = round(XStimParams.loabl:...
        (XStimParams.hiabl - XStimParams.loabl)/(XStimParams.numabls-1):...
        XStimParams.hiabl);
    set(H.ablvector,'string',num2str(XStimParams.abls));
else
    XStimParams.abls = str2num(get(H.ablvector,'string'));    
    XStimParams.numabls = length(XStimParams.abls);
    set(H.lowabl,'String','');
    set(H.highabl,'String','');
    set(H.stepabl,'String','');
end

set(H.numabls,'string',num2str(XStimParams.numabls));
set(H.recorddata_FN, 'String', FN.data);

%Set use loc_flag
if get(H.abl_useHRTF_flag(1),'value')  & XStimParams.loc_flag ~=1    % std hrtfs
    XStimParams.loc_flag = 1;     
    set(H.abl_useHRTF_flag(2),'value',0)
elseif get(H.abl_useHRTF_flag(2),'value') & XStimParams.loc_flag ~= 2     % eq hrtfs
    XStimParams.loc_flag = 2;
    set(H.abl_useHRTF_flag(1),'value',0)
else
    XStimParams.loc_flag =0;
end 

% save current parameter values
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim