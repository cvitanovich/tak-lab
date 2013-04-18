function logplot(H,samplerate,gtitle,color,minfreq,maxfreq);
% LOGPLOT  Plot magnitude of one or more complex signals.
%          MAGPLOT(H,samplerate,gtitle,color,minfreq,maxfreq);
%          required input parameters:
%            H          --> each column of matrix H is plotted
%                           (H should contain samples at frequencies
%                           0 ... samplerate/2)
%            samplerate --> sampling frequency for data
%          optional input parameters:
%            gtitle   --> title of graph (default empty)
%            color    --> the i-th column of H is plotted with the
%                         color given by the i-th entry of this vector
%                         (default yrgbwmc)
%            minfreq  --> lower frequency bound (default 0)
%            maxfreq  --> upper frequency bound (default samplerate/2)
%          The y-axis is labeled "magnitude/dB", the x-axis is labeled
%          "frequency/Hz".

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 20.06.1994 A. Raab

maxcolor = 'yrgbwmc';         % available colors

if nargin<6                   % set parameteres which are not specified
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
while numcolors < size(H,2)                       % make sure that there
  [color,numcolors]=pltstyle([color maxcolor]);   % are enough colors
end;

xtics = linspace(0,samplerate/2,size(H,1));       % tics for x-axis
logtics = log10(xtics/minfreq);

r=log10(maxfreq/minfreq);
%f=minfreq*10.^(linspace(0,r,size(H,1)));
f=minfreq*10.^(linspace(0,r,500));

mag = 20*log10(abs(H(:,1)));                      % first column magnitude
nmag= interp1(xtics,mag,f);
figure;
%semilogx(xtics,mag,setstr(color(1:3)));
semilogx(f,nmag,setstr(color(1:3)));

hold on;                                   % plot next data to same graph
for i=2:size(H,2)
  mag = 20*log10(abs(H(:,i)));
nmag= interp1(xtics,mag,f);
 %semilogx(xtics,mag,setstr(color(3*i-2:3*i)));
h=semilogx(f,nmag,setstr(color(3*i-2:3*i)));
set(h,'LineWidth',2)
end;
hold off;
title(gtitle);
xlabel('frequency/Hz');
ylabel('magnitude/dB');
axisscaling = axis;
axisscaling(1) = minfreq;
axisscaling(2) = maxfreq;
axis(axisscaling);
