function p_layout(file_format,filetype,filename,n_o_subplot,headline,n_o_direction,...
                  direction_matrix,xaxis,n_lines)

% function p_layout(file_format,filetype,filename,n_o_subplot,headline,n_o_direction,...
%                   direction_matrix,xaxis,n_lines)
%
% P_LAYOUT generates a layout for the plots of the phase of transfer functions.
%          Different layouts are provided, depending on the different file formats.
%          The following attributes are controlled:
%
%          - controls the correct range for the x-axis
%          - labels the axes
%          - displayes title
%          - displayes the directions
%
%          meaning of the input arguments:
%          - file_format: must be MTLSIG or SIGMAP, and defines the kind of
%                         labeling
%          - filetype: number of identification for the signal type
%          - filename: name of the file
%          - n_o_subplot: number of subplot in the actual figure (only the
%                         first subplot has a title)
%          - headline: title of the first subplot
%          - n_o_direction: number of direction; index about all plotted
%                           signals
%          - direction_matrix: matrix including the elevation and azimuth
%                              angle of the different directions. The first
%                              row includes the elevation, the second row
%                              the azimuth. Each column represents a
%                              direction.
%          - xaxis: vector including the points of sampling
%          - n_lines: number of lines
%
%          For further information see: ph_plot
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


% axis changing
axissetting = axis;
axissetting(1) = xaxis(1);
axissetting(2) = xaxis(n_lines);
axis(axissetting);

% title and labels
if ( (filetype>=3)&(filetype<=10) )                  % frequency signal
 if n_o_subplot==1
  if ( (filetype==3)|(filetype==4)|(filetype==7)|(filetype==8)|...
       (filetype==9)|(filetype==10) )
    title([headline,'phase  ',filename]);
  elseif (filetype==5)
    title([headline,'phase delay ',filename]);
  elseif (filetype==6)
    title([headline,'group delay ',filename]);
  end;
 end;
 ylabel('phase / deg');
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
 end;
else
 error('choose another layout function');
end

% direction display
if file_format=='MTLSIG'
  dirlabel(filetype,n_o_direction,direction_matrix);
end;