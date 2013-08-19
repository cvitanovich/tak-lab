function zoom_sub(xaxis)

% function zoom_sub(xaxis)
%
% ZOOM_SUB is an user menu for selecting ranges in a figure including more than
%          one subplot.
%          If MATLAB runs under WINDOWS the following menu is attached 
%          to a figure:
%
%          Set X-Range:         use mouse to select displayed x-range
%          Set Y-Range:         use mouse to select displayed y-range
%          Select XY-Range:     use mouse to select displayed x- and y-range
%          Autoscale:           autoscale plot
%          Grid on/off:         set grid on/off
%
%          Before selecting a range, it is necessary to choose the subplot,
%          where the selection should be done. The selected subplot is
%          marked in a different color. The selected range is used for all
%          subplots.
%          The ZOOM_SUB - function is only available for at most four
%          subplots.
%
%          The input argument 'xaxis' is a vector including the points of
%          sampling.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

%if computer == 'PCWIN'

  axes_h = findax(gcf);
  n_o_ax = length(axes_h);
  if n_o_ax>4
    disp('Axis-Control only possible for less than 5 subplots');
  else
    for i=1:n_o_ax
      id = n_o_ax + 1 - i;
      idaxes(axes_h(i),id);
    end;

    item = uimenu('Label','Axis-Control');
    xsub = uimenu(item,'Label','Set X-Range ');
    ysub = uimenu(item,'Label','Set Y-Range ');
    xysub = uimenu(item,'Label','Set XY-Range');
    uimenu(item,'Label','Autoscale','CallBack','zoom_(4,xaxis)');
    uimenu(item,'Label','Grid on/off','Callback','zoom_(5,xaxis)');

    % submenus

    if n_o_ax==1
	uimenu(xsub,'Label','Subplot 1','callback','zoom_(1,xaxis,1)');
	uimenu(ysub,'Label','Subplot 1','callback','zoom_(2,xaxis,1)');
	uimenu(xysub,'Label','Subplot 1','callback','zoom_(3,xaxis,1)');
    elseif n_o_ax==2
	uimenu(xsub,'Label','Subplot 1','callback','zoom_(1,xaxis,1)');
	uimenu(xsub,'Label','Subplot 2','callback','zoom_(1,xaxis,2)');
 
	uimenu(ysub,'Label','Subplot 1','callback','zoom_(2,xaxis,1)');
	uimenu(ysub,'Label','Subplot 2','callback','zoom_(2,xaxis,2)');
  
	uimenu(xysub,'Label','Subplot 1','callback','zoom_(3,xaxis,1)');
	uimenu(xysub,'Label','Subplot 2','callback','zoom_(3,xaxis,2)');
     elseif n_o_ax==3
	uimenu(xsub,'Label','Subplot 1','callback','zoom_(1,xaxis,1)');
	uimenu(xsub,'Label','Subplot 2','callback','zoom_(1,xaxis,2)');
	uimenu(xsub,'Label','Subplot 3','callback','zoom_(1,xaxis,3)');
 
	uimenu(ysub,'Label','Subplot 1','callback','zoom_(2,xaxis,1)');
	uimenu(ysub,'Label','Subplot 2','callback','zoom_(2,xaxis,2)');
	uimenu(ysub,'Label','Subplot 3','callback','zoom_(2,xaxis,3)');

	uimenu(xysub,'Label','Subplot 1','callback','zoom_(3,xaxis,1)');
	uimenu(xysub,'Label','Subplot 2','callback','zoom_(3,xaxis,2)');
	uimenu(xysub,'Label','Subplot 3','callback','zoom_(3,xaxis,3)');
     elseif n_o_ax==4
	uimenu(xsub,'Label','Subplot 1','callback','zoom_(1,xaxis,1)');
	uimenu(xsub,'Label','Subplot 2','callback','zoom_(1,xaxis,2)');
	uimenu(xsub,'Label','Subplot 3','callback','zoom_(1,xaxis,3)');
	uimenu(xsub,'Label','Subplot 4','callback','zoom_(1,xaxis,4)');

	uimenu(ysub,'Label','Subplot 1','callback','zoom_(2,xaxis,1)');
	uimenu(ysub,'Label','Subplot 2','callback','zoom_(2,xaxis,2)');
	uimenu(ysub,'Label','Subplot 3','callback','zoom_(2,xaxis,3)');
	uimenu(ysub,'Label','Subplot 4','callback','zoom_(2,xaxis,4)');

	uimenu(xysub,'Label','Subplot 1','callback','zoom_(3,xaxis,1)');
	uimenu(xysub,'Label','Subplot 2','callback','zoom_(3,xaxis,2)');
	uimenu(xysub,'Label','Subplot 3','callback','zoom_(3,xaxis,3)');
	uimenu(xysub,'Label','Subplot 4','callback','zoom_(3,xaxis,4)');
     end;
  end; % n_o_ax<=4
%end;  % PCWIN

