function [] = SearchCells()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%SearchCells: the GUI for conducting searches for cells
%*******************************************************************************
%	The Graphical User Interface
%*******************************************************************************
XStimParams.curr_stimdur = 100;

F.pos_incr = 20;F.tot_ht = 0;
%Figure window
if(~exist1('H.searchfig') &...
      ~exist1('H.ablfig') &...
      ~exist1('H.itdfig') &...
      ~exist1('H.ildfig') &...
      ~exist1('H.ildfreqfig') &...
      ~exist1('H.spacefig') &...
      ~exist1('H.ildalonefig'))
   Xstim_control_pos = get(H.fig,'Position');
   searchfig_l = 100; searchfig_wd = 350; searchfig_ht = 800;
   H.searchfig = figure('Units','pixels',...
      'Position',[searchfig_l 0 searchfig_wd searchfig_ht],...
      'Name','SEARCH Control',...
      'NumberTitle','off');

F.tot_ht = F.tot_ht + F.pos_incr;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[.25*searchfig_wd searchfig_ht-F.tot_ht 175 F.pos_incr],...
   'ForegroundColor','blue',...
   'FontWeight','bold',...
   'String', 'SEARCH Control');


F.tot_ht = F.tot_ht + F.pos_incr;
% toggle stimulus noise / tone
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[searchfig_wd-140 searchfig_ht-F.tot_ht 60 F.pos_incr],...
   'BackgroundColor',[.8 .8 .8],...
   'String', 'change to:');
H.togglestim = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[searchfig_wd-80 searchfig_ht-F.tot_ht 75 F.pos_incr],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'tone');

F.tot_ht = F.tot_ht + 2*F.pos_incr;
%ITD: slider
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht 80 F.pos_incr],...
   'FontSize',12,...
   'FontWeight','bold',...
   'String', 'ITD');
F.tot_ht = F.tot_ht + F.pos_incr*1.5;
H.ITDslider = uicontrol('Style','slider',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht searchfig_wd-60 F.pos_incr*1.5],...
   'BackgroundColor','blue',...
   'Min',-250,'Max',250,...
   'SliderStep',[5/500 20/500],...
   'Value',XStimParams.curr_ITD,...
   'Callback','setInfo_search');

H.ITDslider_text = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[searchfig_wd-60 searchfig_ht-F.tot_ht 60 F.pos_incr*1.5],...
   'FontSize',12,...
   'FontWeight','bold',...
   'String', [num2str(get(H.ITDslider,'Value')) ' us']);

%ILD: slider
F.tot_ht = F.tot_ht + 2*F.pos_incr;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht 80 F.pos_incr],...
   'FontSize',12,...
   'FontWeight','bold',...
   'String', 'ILD');
F.tot_ht = F.tot_ht + 1.5*F.pos_incr;
H.ILDslider = uicontrol('Style','slider',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht searchfig_wd-60 F.pos_incr*1.5],...
   'BackgroundColor','blue',...
   'Min',-40,'Max',40,...,
   'Value',XStimParams.curr_ILD,...
   'SliderStep',[2/80 10/80],...
   'Callback','setInfo_search');

H.ILDslider_text = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[searchfig_wd-60 searchfig_ht-F.tot_ht 60 F.pos_incr*1.5],...
   'FontSize',12,...
   'FontWeight','bold',...
   'String', [num2str(get(H.ILDslider,'Value')) ' dB']);

%ABL: slider
F.tot_ht = F.tot_ht + 2*F.pos_incr;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht 80 F.pos_incr],...
   'FontWeight','bold',...
   'FontSize',12,...
   'String', 'ABL');
F.tot_ht = F.tot_ht + 1.5*F.pos_incr;
H.ABLslider = uicontrol('Style','slider',...
   'Units','pixels',...
   'BackgroundColor','blue',...
   'Min',-119,'Max',-20,...
   'Value',XStimParams.curr_ABL,...
   'SliderStep',[2/100 10/100],...
   'Position',[0 searchfig_ht-F.tot_ht searchfig_wd-60 F.pos_incr*1.5],...
   'Callback','setInfo_search');
H.ABLslider_text = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[searchfig_wd-60 searchfig_ht-F.tot_ht 60 F.pos_incr*1.5],...
   'FontSize',12,...
   'FontWeight','bold',...
   'String', [num2str(get(H.ABLslider,'Value')) ' dB']);


F.tot_ht = F.tot_ht + F.pos_incr/2;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht+F.pos_incr/2 searchfig_wd-10 F.pos_incr/2],...
   'FontSize',12,...
   'String', '|               |                |                 |                |                 |               ');
F.tot_ht = F.tot_ht + F.pos_incr/2;
H.SPLtxt = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht searchfig_wd-10 F.pos_incr],...
   'FontSize',12,...
   'String', '-59          -42           -20             0             18             38          dB SPL');
 
%ISI: edit box
F.tot_ht = F.tot_ht + 3*F.pos_incr;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht 80 F.pos_incr],...
   'String', 'ISI (ms)');
H.ISI = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[90 searchfig_ht-F.tot_ht 100 F.pos_incr],...
   'BackgroundColor',[.7 .7 .7],...
   'Value',XStimParams.search_ISI,...
   'string',num2str(XStimParams.search_ISI),...
   'Callback','setInfo_search');

%DUR: edit box
F.tot_ht = F.tot_ht + 2*F.pos_incr;
uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht 80 F.pos_incr],...
   'String', 'DUR (ms)');
H.DUR = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[80 searchfig_ht-F.tot_ht 100 F.pos_incr],...
   'BackgroundColor',[.7 .7 .7],...
   'Value',XStimParams.curr_stimdur,...
   'string',num2str(XStimParams.curr_stimdur),...
   'Callback','setInfo_search');

%Frequency for tones
F.tot_ht = F.tot_ht + 3*F.pos_incr;
H.tonal_frequency_title = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht 100 F.pos_incr],...
   'FontWeight','bold',...
   'FontSize',12,...
   'String', 'Frequency');

F.tot_ht = F.tot_ht + 1.5*F.pos_incr;
H.tonal_frequency = uicontrol('Style','slider',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht searchfig_wd-110 F.pos_incr*1.5],...
   'Min',2000,...
   'Max',11000,...
   'SliderStep',[100/9000 1000/9000],...
   'Value',XStimParams.curr_freq,...
   'Callback','setInfo_search');
H.tonal_frequency_text = uicontrol('Style','edit',...
   'Units','pixels',...
   'Position',[searchfig_wd-95 searchfig_ht-F.tot_ht 55 F.pos_incr*1.5],...
   'String',num2str(get(H.tonal_frequency,'Value')),...
   'FontSize',12,...
   'FontWeight','bold',...
   'Callback','setInfo_search');
H.tonal_frequency_hz = uicontrol('Style','text',...
   'Units','pixels',...
   'Position',[searchfig_wd-40 searchfig_ht-F.tot_ht 30 F.pos_incr],...
   'FontSize',12,...
   'FontWeight','bold',...
   'String','Hz');


% modulate sound
    F.tot_ht = F.tot_ht + 2* F.pos_incr;
    switch XStimParams.mod_type
        case 'Tone'
            mod_num = 1;
        case 'Sq wave'
            mod_num = 2;
        case 'LP Noise'
            mod_num = 3;
        case 'File'
            mod_num = 4;
        case 'None'
            mod_num = 5;
        otherwise
            mod_num = 5;
    end
    uicontrol('Style','text',...
        'Units','pixels',...
        'FontWeight','bold',...
        'Position',[0 searchfig_ht-F.tot_ht 70 F.pos_incr],...
        'String', ' AM? ');
    H.Search_mod_type = uicontrol('Style','popup',...
        'Units','pixels',...
        'Position',[70 searchfig_ht-F.tot_ht 70 F.pos_incr],...
        'BackgroundColor','White',...
        'String',...
        'Tone|Sq wave|LP Noise|File|None',...
        'Value',mod_num,...
        'Callback','setInfo_search');
    H.Search_mod_txt = uicontrol('Style','text',...
        'Units','pixels',...
        'HorizontalAlignment','right',...
        'Position',[140 searchfig_ht-F.tot_ht 40 F.pos_incr],...
        'String', 'freq  ');
    H.Search_mod_freq = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[180 searchfig_ht-F.tot_ht 40 F.pos_incr],...
        'String', num2str(XStimParams.mod_freq(1)),...
        'Callback','setInfo_search');
    
    H.Search_mod_txtA = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[220 searchfig_ht-F.tot_ht 30 F.pos_incr],...
        'String', 'depth');
    H.Search_mod_depth = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[250 searchfig_ht-F.tot_ht 40 F.pos_incr],...
        'String', num2str(XStimParams.mod_depth(1)),...
        'Callback','setInfo_search');
    H.Search_mod_txtB = uicontrol('Style','text',...
        'Units','pixels',...
        'Position',[290 searchfig_ht-F.tot_ht 30 F.pos_incr],...
        'String', 'phase');
    H.Search_mod_phase = uicontrol('Style','edit',...
        'Units','pixels',...
        'Position',[320 searchfig_ht-F.tot_ht 30 F.pos_incr],...
        'String', num2str(XStimParams.mod_phase(1)),...
        'Callback','setInfo_search');
    
    H.Search_mod_pb = uicontrol('Style','pushbutton',...
        'Units','pixels',...
        'Position',[140 searchfig_ht-F.tot_ht 25 F.pos_incr],...
        'BackgroundColor','red',...
        'ForegroundColor','blue',...
        'String', 'FN',...
        'Visible','off',...
        'Callback','setInfo_search');

%Plotting
F.tot_ht = F.tot_ht + 6* F.pos_incr;
H.plotpsd = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[searchfig_wd-100 searchfig_ht-F.tot_ht 95 F.pos_incr],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Plot PSD');

%Pause
H.pausesearch = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht 80 F.pos_incr],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'Pause');


% toggle to clear plot
H.clearplot_toggle = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[searchfig_wd-200 searchfig_ht-F.tot_ht 75 F.pos_incr],...
   'BackgroundColor','blue',...
   'ForegroundColor','white',...
   'FontWeight','bold',...
   'String', 'clearPlot');


%Exit
F.tot_ht = F.tot_ht + 2* F.pos_incr;
H.exitsearch = uicontrol('Style','togglebutton',...
   'Units','pixels',...
   'Position',[0 searchfig_ht-F.tot_ht 80 F.pos_incr],...
   'BackgroundColor','red',...
   'ForegroundColor','yellow',...
   'FontWeight','bold',...
   'String', 'EXIT',...
   'Callback','Exit_Search');

end %end GUI

setInfo_search;

engage_search