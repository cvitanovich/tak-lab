function zoom_(id,xaxis,n_o_axes)

% function zoom_(id,xaxis,n_o_axes)
%
% ZOOM_ is called from ZOOM_SUB.
%       Parameter 'id' selects the desired action:
%       id   action             explanation
%       01   Set X-Range        use mouse to select displayed x-range
%       02   Set Y-Range        use mouse to select displayed y-range
%       03   Select Window      use mouse to select displayed x- and y-range
%       04   Autoscale          autoscale plot
%       05   Grid on/off        set grid on/off
%
%       Parameter 'xaxis' is a vector including the points of sampling.
%       Parameter 'n_o_axes' includes the number of subplots in the actual
%       figure. This zoom function is available in figures that may include
%       more than one subplot.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


axes_h = findax(gcf);
if (id<=3)
  for i=1:length(axes_h)
    p=get(get(axes_h(i),'title'),'Userdata');
    if p==n_o_axes
     handle_nr = i;
     break;
    end;
  end;
  set(gcf,'CurrentAxes',axes_h(handle_nr));
  set(gca,'XColor','y','YColor','y');
end;

if id == 1
  axissetting = axis;           % get current axis setting
  [x1,y1] = ginput(1);          % get lower x-value
  [x2,y2] = ginput(1);          % get upper x-value
  if x1 > x2                    % swap if necessary
    y1=x1; x1=x2; x2=y1;
  end;
  axissetting(1) = x1;          % set new x-values
  axissetting(2) = x2;
  set(gca,'XColor','w','YColor','w');
  for i=1:length(axes_h);                     % change all figure axes elements
    axes(axes_h(i));
    axis(axissetting);
  end
elseif id == 2
  axissetting = axis;
  [x1,y1] = ginput(1);
  [x2,y2] = ginput(1);
  if y1 > y2
    x1=y1; y1=y2; y2=x1;
  end;
  axissetting(3) = y1;
  axissetting(4) = y2;
  set(gca,'XColor','w','YColor','w');
  for i=1:length(axes_h)                     % change all figure axes elements
    axes(axes_h(i));
    axis(axissetting);
  end;
elseif id == 3
  axissetting = axis;
  [x1,y1] = ginput(1);
  [x2,y2] = ginput(1);
  if x1>x2 d=x1; x1=x2; x2=d; end;
  if y1>y2 d=y1; y1=y2; y2=d; end;
  set(gca,'XColor','w','YColor','w');
  for i=1:length(axes_h);                     % change all figure axes elements
    axes(axes_h(i));
    axis([x1 x2 y1 y2]);
  end;
elseif id == 4
  axis('auto');
  axissetting = axis;
  axissetting(1) = xaxis(1);
  axissetting(2) = xaxis(length(xaxis));
  for i=1:length(axes_h);                     % change all figure axes elements
   axes(axes_h(i));
   axis(axissetting);
  end;
elseif id == 5
  for i=1:length(axes_h);                     % change all figure axes elements
   axes(axes_h(i));
   grid;
  end  ;
end;  
