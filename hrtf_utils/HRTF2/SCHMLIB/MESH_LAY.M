% MESH_LAY generates a layout for 'mesh plots'. Different layouts are provided
%          for the different file formats.
%          The following attributes are controlled:
%
%          - controls the correct range for the x- and y-axis
%          - labels the axes
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

% axis changing

axissetting = axis;
axissetting(1) = xaxis(1);
axissetting(2) = xaxis(n_lines);
axissetting(3) = yaxis(1);
axissetting(4) = yaxis(n_directions);
axis(axissetting);

if filetype==2
  zlabel('h(t)');
elseif z_status==0
  zlabel('A / dB');
else
  zlabel('phase / deg');
end;

if file_format=='MTLSIG'
  if ang==1    
    ylabel('elev,azi / deg');
  else
    ylabel('azi,elev / deg');
  end;
  set(gca,'ytick',yaxis,'yticklabels',lab);
elseif file_format=='SIGMAP'
  ylabel('signal index');
  set(gca,'ytick',yaxis);
end;

if ( (filetype==3)|(filetype==4)|(filetype==5)|(filetype==6))
  xlabel('f / Hz');
elseif (filetype==10)
  xlabel('log10(f) / Hz');
elseif (filetype==7)
  xlabel('f / mel');
elseif (filetype==8)
  xlabel('f / bark');
elseif (filetype==9)
  xlabel('f / ERB');
end;