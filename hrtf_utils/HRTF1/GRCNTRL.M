function grcntrl()
% GRCNTRL User menu for selecting plot windows.
%         If MATLAB runs under WINDOWS the following menu is attached
%         to a plot window:
%         Set X-Range     : use mouse to select displayed x-range
%         Set Y-Range     : use mouse to select displayed y-range
%         Select Window   : use mouse to select displayed x- and y-range
%         Autoscale       : autoscale plot
%         Grid on/off     : set grid on/off

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 20.06.1994 A. Raab

%if computer == 'PCWIN'
  item = uimenu('Label','Axis-Control');
         uimenu(item,'Label','Set X-Range ','Callback','grcntrl_(1)');
         uimenu(item,'Label','Set Y-Range ','Callback','grcntrl_(2)');
         uimenu(item,'Label','Set XY-Range','Callback','grcntrl_(3)');
         uimenu(item,'Label','Autoscale','Callback','grcntrl_(4)');
         uimenu(item,'Label','Grid on/off','Callback','grcntrl_(5)');
%end;
