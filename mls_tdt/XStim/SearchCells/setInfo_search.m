% setInfo_search

%ABL Slider
set(H.ABLslider_text,'String',[num2str(round(get(H.ABLslider,'Value'))) ' dB']);
XStimParams.curr_ABL = round(get(H.ABLslider,'Value'));
set(H.ABLslider,'Value',XStimParams.curr_ABL);

%DUR edit box
XStimParams.curr_stimdur = str2num(get(H.DUR,'string'));
set(H.DUR,'value',XStimParams.curr_stimdur);

%ILD Slider
set(H.ILDslider_text,'String',[num2str(round(get(H.ILDslider,'Value'))) ' dB']);
XStimParams.curr_ILD = get(H.ILDslider,'Value');

%ISI edit box
XStimParams.search_ISI = str2num(get(H.ISI,'string'));
set(H.ISI,'Value',XStimParams.search_ISI);

%ITD Slider
set(H.ITDslider_text,'String',[num2str(round(get(H.ITDslider,'Value'))) ' us']);
XStimParams.curr_ITD = get(H.ITDslider,'Value');

% Tonal Frequency Slider
set(H.tonal_frequency_text,'String',num2str(get(H.tonal_frequency,'Value')));

% Tonal Frequency edit box
textval = str2num(get(H.tonal_frequency_text,'String'));
minval = get(H.tonal_frequency,'Min');
maxval = get(H.tonal_frequency,'Max');
if(textval < minval) textval = minval; end
if(textval > maxval) textval = maxval; end
set(H.tonal_frequency,'Value',textval);
set(H.tonal_frequency_text,'String',num2str(textval));
XStimParams.curr_freq = get(H.tonal_frequency,'Value');

% check out modulation parameters
XStimParams.mod_depth(1) = str2num(get(H.Search_mod_depth,'String'));
mod_type = get(H.Search_mod_type,'String');
mod_num = get(H.Search_mod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Tone'
        set(H.Search_mod_txt,'String', 'freq  ');
        set(H.Search_mod_freq,'Visible','on');
        set(H.Search_mod_phase,'Visible','on');
        set(H.Search_mod_txtA,'Visible','on');
        set(H.Search_mod_depth,'Visible','on');
        set(H.Search_mod_txtB,'Visible','on');
        set(H.Search_mod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Search_mod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.Search_mod_phase,'String'));
    case 'File'
        set(H.Search_mod_pb,'Visible','on');
        set(H.Search_mod_freq,'Visible','off');
        set(H.Search_mod_txtA,'Visible','on');
        set(H.Search_mod_depth,'Visible','on');
        set(H.Search_mod_phase,'Visible','off');
        set(H.Search_mod_txtB,'Visible','off');
        if get(H.Search_mod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.Search_mod_pb,'Value',0);
        end
        set(H.Search_mod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.Search_mod_txt,'String', 'CutOff ');
        set(H.Search_mod_freq,'Visible','on');
        set(H.Search_mod_txtA,'Visible','on');
        set(H.Search_mod_depth,'Visible','on');
        set(H.Search_mod_phase,'Visible','off');
        set(H.Search_mod_pb,'Visible','off');
        set(H.Search_mod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.Search_mod_freq,'String'));
        
    case 'None'
        set(H.Search_mod_txt,'String', 'no mod  ');
        set(H.Search_mod_freq,'Visible','off');
        set(H.Search_mod_phase,'Visible','off');
        set(H.Search_mod_pb,'Visible','off');
        set(H.Search_mod_txtB,'Visible','off');
        set(H.Search_mod_txtA,'Visible','off');
        set(H.Search_mod_depth,'Visible','off');
        
    otherwise
end

