function phasplot(H,samplerate,gtitle,color,minfreq,maxfreq,correct,unwraptol);
% PHASPLOT Plot phase of one or more complex signals in degrees.
%          PHASPLOT(H,samplerate,gtitle,color,minfreq,maxfreq,...
%                   correct,unwraptol);
%          required input parameters:
%            H          --> each column of matrix H is plotted
%                           (H should contain samples at frequencies
%                           0 ... samplerate/2)
%            samplerate --> sampling frequency for data
%          optional input parameters:
%            gtitle    --> title of graph (default empty)
%            color     --> the i-th column of H is plotted with the
%                          color given by the i-th entry of this vector
%                          (default yrgbwmc)
%            minfreq   --> lower frequency bound (default 0)
%            maxfreq   --> upper frequency bound (default samplerate/2)
%            correct   --> phase correction with respect to
%                          the phase of line correct
%            unwraptol --> tolerance for unwrapping (see UNWRAP)
%          The y-axis is labeled "phase/deg", the x-axis is labeled
%          "frequency/Hz".

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 20.06.1994 A. Raab

maxcolor = 'yrgbwmc';      % available colors
docorrect=1;

if nargin<8                % set parameters which are not specified
  unwraptol = pi;
  if nargin<7
    docorrect=0;
    if nargin<6
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
  end;
end;

if ~checkfft(H(:,1))
  disp('warning: H may contain samples over the nyquist frequency');
end;

[color,numcolors]=pltstyle(color);
while numcolors < size(H,2)                     % make sure that there
  [color,numcolors]=pltstyle([color maxcolor]); % are enough colors
end;

xtics = linspace(0,samplerate/2,size(H,1));   % tics for x-axis

phase = 180*unwrap(angle(H(:,1)))/pi;         % first column phase

if docorrect==1                          % phase correction with respect
  pc=phase(correct);                     % to the phase of the first
phase=phase-pc*ones(size(phase))+180*angle(H(correct,1))/pi;

end;
                                         % column of H


figure;
semilogx(xtics,phase,setstr(color(1:3)));
hold on;                                 % plot next data to same graph
for i=2:size(H,2)
  phase = 180*unwrap(angle(H(:,i)))/pi;
  if docorrect==1
  pc=phase(correct);                     % to the phase of the first
  phase=phase-pc*ones(size(phase))+180*angle(H(correct,i))/pi;

%    cv=360*fix((phase(correct)-pc)/360);
%    phase=phase-cv;
  end;
  semilogx(xtics,phase,setstr(color(3*i-2:3*i)));
end;
hold off;
title(gtitle);
xlabel('frequency/Hz');
ylabel('phase/deg');
axisscaling = axis;
axisscaling(1) = minfreq;
axisscaling(2) = maxfreq;
axis(axisscaling);

grcntrl