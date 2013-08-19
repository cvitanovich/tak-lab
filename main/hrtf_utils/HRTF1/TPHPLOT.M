function tphplot(H,samplerate,gtitle,color,minfreq,maxfreq,correct,unwraptol);
% PRELIMINARY
% TPHPLOT Plot phase delay of one or more complex signals in degrees.
%          TPHPLOT(H,samplerate,gtitle,color,minfreq,maxfreq,correct,unwraptol);
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
%            correct    --> sample for phase correction with respect to
%                           the phase of H(1,:)
%            unwraptol  --> tolerance for unwrapping (see UNWRAP)
%          The y-axis is labeled "phase/deg", the x-axis is labeled 
%          "frequency/Hz".

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 20.06.1994 A. Raab 

maxcolor = 'yrgbwmc';           % available colors
docorrect=1;

if nargin<8                     % set parameteres which are not specified
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
while numcolors < size(H,2)
  [color,numcolors]=pltstyle([color maxcolor]); % make sure that thera are enough colors
end;

xtics = linspace(0,samplerate/2,size(H,1));     % tics for x-axis
omega= xtics*2*pi;

phase = unwrap(angle(H(:,1)));          % first column phase
tph=-phase(:)./omega(:);

if docorrect==1                                 % phase correction with respect
  pc=phase(correct);                            % to the phase of the first
end;                                            % column of H

figure;
semilogx(xtics,tph,setstr(color(1:3)));
hold on;                                        % plot next data to same graph
for i=2:size(H,2)
  phase = unwrap(angle(H(:,i)));                
  tph=-phase(:)./omega(:);
  if docorrect==1
    cv=360*fix((phase(correct)-pc)/360);
    phase=phase-cv;
  end;
  semilogx(xtics,tph,setstr(color(3*i-2:3*i)));
end;
hold off;
title(gtitle);
xlabel('frequency/Hz');
ylabel('phase delay/s');
axisscaling = axis;
axisscaling(1) = minfreq;
axisscaling(2) = maxfreq;
axis(axisscaling);
