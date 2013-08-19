% PSC_ING Pseudocolor and contour plot of the actual signal matrix, that is
%         resampled with a lower sampling rate for smoothing the contour
%         lines. If a contour plot of the magnitude is done, the contour
%         lines are selected by choosing the level difference. Contour plots
%         of the phase are done with 10 different level lines.
%
%         Further the following attributes are provided:
%          - controls the correct range for the x-axis and y-axis
%          - labels the axes
%          - automatical labeling of the contour lines
% 
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

resample_factor = 10;
[u_signal,new_xaxis,newlines] = u_samples(signal,xaxis,n_directions,...
                                        n_lines,resample_factor);

if status==0               % contour of the magnitude
  minlevel = fix(min(min(u_signal)));
  maxlevel = fix(max(max(u_signal)));
  disp(['signal mimimum level: ',num2str(minlevel),' dB']);
  disp(['signal maximum level: ',num2str(maxlevel),' dB']);
  dbsteps = [];
  while isempty(dbsteps)      % selecting contour line difference
    dbsteps = input('contour line difference in [db] ? ');
    if dbsteps<0
      dbsteps = -dbsteps;
    end; 
  end;
  minline = dbsteps * fix(minlevel/dbsteps);
  maxline = dbsteps * fix(maxlevel/dbsteps);
  number_of_contours = (maxline-minline)/dbsteps + 1;
  level = linspace(minline,maxline,number_of_contours);
elseif status==1           % contour of the phase
  level = 10;
end;

fig_h = figure;
if file_format=='MTLSIG'
  [yaxis,lab,ang] = y_axis(direction_matrix,n_directions);
elseif file_format=='SIGMAP'
  yaxis = linspace(1,n_directions,n_directions);
end;

pcolor(new_xaxis,yaxis,u_signal);
shading interp;
hold on;
C = contour(new_xaxis,yaxis,u_signal,level);

if scaling=='log'
 set(gca,'xscale','log');
end;

% automatical labeling of the contour lines
if status==0   % contour of the magnitude
 clabel(C,level);
elseif status==1   %contour of the phase
 clabel(C);
end;

% axis changing
axissetting = axis;
axissetting(1) = new_xaxis(1);
axissetting(2) = new_xaxis(newlines);
axissetting(3) = yaxis(1);
axissetting(4) = yaxis(n_directions);
axis(axissetting);

if file_format=='MTLSIG'
  if ang==1    
    ylabel('elev,azi / deg');
  else
    ylabel('azi,elev / deg');
  end;
  set(gca,'TickDir','out','ytick',yaxis,'yticklabels',lab);
elseif file_format=='SIGMAP'
 ylabel('signal index');
 set(gca,'TickDir','out','ytick',yaxis);
end;

if ( (filetype==3)|(filetype==4)|(filetype==5)|(filetype==6))
  xlabel('f / Hz');
elseif (filetype==10)
  xlabel('log10(f / Hz)');
elseif (filetype==7)
  xlabel('f / mel');
elseif (filetype==8)
  xlabel('f / bark');
elseif (filetype==9)
  xlabel('f / ERB');
end;

