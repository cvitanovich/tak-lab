function tgrpplot(H,samplerate,gtitle,color,minfreq,maxfreq);
% TGRPPLOT Plot group delay of one or more complex signals.
%          TGRPPLOT(H,samplerate,gtitle,color,minfreq,maxfreq);
%          required input parameters:
%            H          --> each column of matrix H is plotted
%                           (H should contain samples at frequencies 
%                            0 ... samplerate/2)
%            samplerate --> sampling frequency for data
%          optional input parameters:
%            gtitle     --> title of graph (default empty)
%            color      --> the i-th column of H is plotted with the 
%                           color given by the i-th entry of this vector
%                           (default yrgbwmc)
%            minfreq    --> lower frequency bound (default 0)
%            maxfreq    --> upper frequency bound (default samplerate/2)
%          The y-axis is labeled "group delay/sec", the x-axis is labeled 
%          "frequency/Hz".

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 20.06.1994 A. Raab 

maxcolor = 'yrgbwmc';           % available colors

if nargin<6                     % set parameteres which are not specified
  maxfreq = samplerate/2;       
  if nargin<5                   
    minfreq = 0;  
    if nargin<4
      color=maxcolor;
      if nargin<3
	gtitle='';
      end;
    end;
  end;
end;

if ~checkfft(H(:,1))
  disp('warning: H may contain samples over the nyquist frequency');
end;

[color,numcolors]=pltstyle(color);
while numcolors < size(H,2)
  [color,numcolors]=pltstyle([color maxcolor]); % make sure that thera are enough colors
end;

xtics = linspace(0,samplerate/2,size(H,1));     % tics for x-axis

figure;
semilogx(xtics,tgroup(H(:,1),samplerate),setstr(color(1:3)));
hold on;                                        % plot next data to same graph
for i=2:size(H,2)
  semilogx(xtics,tgroup(H(:,i),samplerate),setstr(color(3*i-2:3*i)));
end;
hold off;
title(gtitle);
xlabel('frequency/Hz');
ylabel('group delay/sec');
axisscaling = axis;
axisscaling(1) = minfreq;
axisscaling(2) = maxfreq;
axis(axisscaling);

