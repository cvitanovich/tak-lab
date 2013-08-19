function v3dcntrl(gf)

% function v3dcntrl(gf)
%
% V3DCNTRL changes the point of view of three-dimensional plots.
%          V3DCNTRL(gf) attaches itself to the 3D axis in the figure
%          window specified by the handle 'gf'. If 'gf' is omitted, the
%          current figure is used.
%
%          V3DCNTRL offers the possibility to find the best orientation
%          for a three dimensional plot. It can be used in figure windows
%          which have a single 3D axes. Upon its use, V3DCNTRL repositions
%          the axes, leaving room for the following controls:
%
%           - A scroll bar controlling the azimuth.
%           - A scroll bar controlling the elevation.
%           - "Azimuth" and "Elevation" edit boxes to display
%             the current values of the associated scroll bars
%             or allow the desired azimuth or elevation to be
%             entered directly.
%
%           Notes: Only one V3DCNTRL per figure window is allowed.
%                  Only one axis can be present in the figure window.
%
%           The function V3DCNTRL is based on the function V3DTOOL by
%           Denis W. Brown in the MATLAB zoomtool toolbox.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

% check and see if progtools are installed

%if computer=='PCWIN'

 if exist('findpush') ~= 2,
    error('v3dtool: Requires D. W. Brown''s programming toolbox...');
 end;
 % check input
 if nargin > 1,
    error('v3dtool: Invalid number of input arguments...');
 end;

 % get current axis if none given
 if nargin ~= 1,
	gf = gcf;
 end;
 % check for only one axes
 if cntaxes(gf) ~= 1,
    error('v3dtool: Only one axes per figure allowed...');
 end;
 % check if already running in window
 if finduitx(gf,'Elevate'),
    error('v3dtool: Already running in figure window...');
 end;
 % get handle figure axes
 ga = gca;

 % check axis
 if isempty(get(ga,'ZTickLabels')),
    error('v3dtool: Axes is not a 3D axes...');
 end;

 % don't allow any other plots until quit is pushed (so the 
 %   controls are removed from the window)
 set(gf,'NextPlot','new');

 % set standard azimuth and elevation
 v = get(ga,'View');
 if v(1)==322.5;
   v(1) =  -37.5;
 end;

 % get the axes limits and save it away
 p = get(ga,'Position');
 set(ga,'UserData',p);

 % change axis size to fit controls in window
 p = [.15 .15 .7 .75];
 set(ga,'position',p);

 % redraw command
 redraw = 'view(az,el);';

 % get slider values
 getslider = [...
    'az = get(findslid(gcf,''azimuth''),''Value'');'...
    'el = get(findslid(gcf,''elevation''),''Value'');'...
	     ];

 % get slid edit box values
 getsledit = [...
    'az = str2num(get(findedit(gcf,''azimuth''),''String''));'...
    'el = str2num(get(findedit(gcf,''elevation''),''String''));'...
             ];

 % set slider values
 setslider = [...
	'set(findslid(gcf,''azimuth''),''Value'',az);'...
	'set(findslid(gcf,''elevation''),''Value'',el);'...
	     ];
 % set slider edit boxes
 setsledit = [...
	'set(findedit(gcf,''azimuth''),''String'',num2str(az));'...
	'set(findedit(gcf,''elevation''),''String'',num2str(el));'...
             ];

 % check edit box values, if in error, reset to present values,
 %   do not redraw, otherwise redraw
 chkedit = [...
    'if (abs(az) > 180) | (el < -90) | (el > 90),'...
        getslider ...
        setsledit ...
    'else,' ...
        getsledit ...
        setslider ...
        redraw ...
    'end;'...
           ];

 % slider callback
 slidecall = [getslider setsledit redraw];

 % edit callback
 sleditcall = [getsledit chkedit];

 %slider position
 sposa = [.1 0 .8 .04];      % azimuth
 spose = [.95 .1 .04 .8];    % elevation

 % make slider contols
 h = uicontrol(gf,'Style','slider','Units','Normal',...
	'Position',sposa,'String','azimuth',...
	'UserData',ga,'Max',180,'Min',-180,'Value',v(1),...
	'Callback',slidecall);
 set(h,'Value',v(1));
 h =uicontrol(gf,'Style','Slider','Units','Normal',...
	'Position',spose,'String','elevation',...
	'UserData',ga,'Max',90,'Min',-90,...
	'Callback',slidecall);
 set(h,'Value',v(2));
 % label sliders
 cback = get(gf,'Color');
 cfore = get(ga,'XColor');
 uicontrol(gf,'Style','Text','Units','Normal',...
    'Position',[sposa(1)-0.06 sposa(2) 0.06 sposa(4)]',...
    'BackGround',cback,'ForeGround',cfore,...
    'String','-180','Horiz','left');
 uicontrol(gf,'Style','Text','Units','Normal',...
    'Position',[sposa(1)+sposa(3)+0.005 sposa(2) 0.07 sposa(4)]',...
    'BackGround',cback,'ForeGround',cfore,...
    'String','+180','Horiz','left');
 uicontrol(gf,'Style','Text','Units','Normal',...
    'Position',[spose(1) spose(2)-0.05 spose(3) 0.05]',...
    'BackGround',cback,'ForeGround',cfore,...
    'String','-90','Horiz','center');
 uicontrol(gf,'Style','Text','Units','Normal',...
    'Position',[spose(1)-0.01 spose(2)+spose(4) spose(3)+0.04 0.05]',...
    'BackGround',cback,'ForeGround',cfore,...
    'String','+90','Horiz','center');

 % edit boxes for setting azimuth and elevation
 cback = 'blue';
 uicontrol(gf,'Style','edit','Units','normal','UserData','azimuth',...
    'Position',[sposa(1)+0.15 sposa(2)+sposa(4)+0.01 0.1 0.04],...
    'BackGround',cback,'ForeGround',cfore,...
    'String',num2str(v(1)),'Horiz','right',...
    'CallBack',sleditcall);
 uicontrol(gf,'Style','edit','Units','normal','UserData','elevation',...
    'Position',[spose(1)-0.15 spose(2)+spose(4)-0.05 0.13 0.04],...
    'BackGround',cback,'ForeGround',cfore,...
    'String',num2str(v(2)),'Horiz','right',...
    'CallBack',sleditcall);

 % labels for text objects
 cback = get(gf,'Color');
 p = get(findedit(gf,'azimuth'),'Position');
 uicontrol(gf,'Style','text','Units','normal','UserData','azimuth',...
    'Position',[p(1)-0.15 p(2) 0.12 p(4)],...
    'BackGround',cback,'ForeGround',cfore,...
    'String','Azimuth:','Horiz','center');
 p = get(findedit(gf,'elevation'),'Position');
 uicontrol(gf,'Style','text','Units','normal','UserData','elevation',...
    'Position',[p(1) p(2)+p(4) p(3) p(4)],...
    'BackGround',cback,'ForeGround',cfore,...
    'String','Elevation','Horiz','center');

end;
