function zoom1_(id,xaxis)

% function zoom1_(id,xaxis)
%
% ZOOM1_ Is called from ZOOM1. 
%        Parameter id selects the desired action:
%        id   action            explanation
%        01   Set X-Range       use mouse to select displayed x-range
%        02   Set Y-Range       use mouse to select displayed y-range
%        03   Set XY-Range      use mouse to select displayed x- and y-range
%        04   Autoscale         autoscale plot
%        05   Grid on/off       set grid on/off
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


if id == 1
  axissetting = axis;           % get current axis setting
  [x1,y1] = ginput(1);          % get lower x-value
  [x2,y2] = ginput(1);          % get upper x-value
  if x1 > x2                    % swap if necessary
    y1=x1; x1=x2; x2=y1;
  end;
  axissetting(1) = x1;          % set new x-values
  axissetting(2) = x2;
  axis(axissetting);
elseif id == 2
  axissetting = axis;
  [x1,y1] = ginput(1);
  [x2,y2] = ginput(1);
  if y1 > y2
    x1=y1; y1=y2; y2=x1;
  end;
  axissetting(3) = y1;
  axissetting(4) = y2;
  axis(axissetting);
elseif id == 3
  axissetting = axis;
  [x1,y1] = ginput(1);
  [x2,y2] = ginput(1);
  if x1>x2 d=x1; x1=x2; x2=d; end;
  if y1>y2 d=y1; y1=y2; y2=d; end;
  axis([x1 x2 y1 y2]);
elseif id == 4
  axis('auto');
  axissetting = axis;
  axissetting(1) = xaxis(1);
  axissetting(2) = xaxis(length(xaxis));
  axis(axissetting);
elseif id == 5
  grid;
end;  
