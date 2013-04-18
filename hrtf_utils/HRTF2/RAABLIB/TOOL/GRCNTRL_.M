function grcntrl_(id)
% GRCNTRL_ Is called from GRCNTRL.
%          Parameter id selects the desired action:
%          id   action          explanation
%          01   Set X-Range     use mouse to select displayed x-range
%          02   Set Y-Range     use mouse to select displayed y-range
%          03   Select Window   use mouse to select displayed x- and y-range
%          04   Autoscale       autoscale plot
%          05   Grid on/off     set grid on/off
%          06   display coordinates 

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 20.06.1994 A. Raab

global X NUMBER_OF_LINES

if id == 1
  axissetting = axis;          % get current axis setting
  [x1,y1] = ginput(1);         % get lower x-value
  [x2,y2] = ginput(1);         % get upper x-value
  if x1 > x2                   % swap if necessary
    y1=x1; x1=x2; x2=y1;
  end;
  axissetting(1) = x1;         % set new x-values
  axissetting(2) = x2;
  axis(axissetting);
end
else
if id == 2
  axissetting = axis;
  [x1,y1] = ginput(1);
  [x2,y2] = ginput(1);
  if y1 > y2
    x1=y1; y1=y2; y2=x1;
  end;
  axissetting(3) = y1;
  axissetting(4) = y2;
  axis(axissetting);
end
else
if id == 3
  axissetting = axis;
  [x1,y1] = ginput(1);
  [x2,y2] = ginput(1);
  if x1>x2 d=x1; x1=x2; x2=d; end;
  if y1>y2 d=y1; y1=y2; y2=d; end;
  axis([x1 x2 y1 y2]);
end
else
if id == 4
  axis('auto');
end
else
if id == 5
 grid;
end;
else 
if id==6
[x,y]=ginput(1);
disp(['X= ' num2str(x) ' Hz     Y= ' num2str(y) ' dB']);  
end;
