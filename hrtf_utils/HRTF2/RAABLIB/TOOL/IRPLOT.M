function irplot(h,samplerate,gtitle,color,mintime,maxtime);
% IRPLOT   Plot one or more real (time) signals.
%          IRPLOT(h,samplerate,gtitle,color,mintime,maxtime);
%          required input parameters:
%            h          --> each column of matrix h is plotted
%            samplerate --> sampling frequency for data
%          optional input parameters:
%            gtitle     --> title of graph (default empty)
%            color      --> the i-th column of h is plotted with the
%                           color given by the i-th entry of this vector
%                           (default yrgbwmc)
%            mintime    --> lower time limit in milliseconds (default 0)
%            maxfreq    --> upper time limit in milliseconds (default all)
%          The y-axis is labeled "h(t)", the x-axis is labeled
%          "time/ms".

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.06.1994 A. Raab

T=1000/samplerate;              % sample time in milliseconds
maxcolor = 'yrgbwmc';           % available colors

if nargin<6                     % set parameteres which are not specified
  maxtime = T*(size(h,1)-1);
  if nargin<5
    mintime = 0;
    if nargin<4
      color=maxcolor;
      if nargin<3
        gtitle='';
      end;
    end;
  end;
end;

while length(color) < size(h,2)
  color = [color maxcolor];     % make sure that thera are enough colors
end;

xtics=0:T:(T*(size(h,1)-1));

figure;
hold on;                        % plot all data to same graph
for i=1:size(h,2)
  plot(xtics,h(:,i),setstr(color(i)));
end;

title(gtitle);
xlabel('time/ms');
ylabel('h(t)');
axisscaling = axis;
axisscaling(1) = mintime;
axisscaling(2) = maxtime;
axis(axisscaling);
