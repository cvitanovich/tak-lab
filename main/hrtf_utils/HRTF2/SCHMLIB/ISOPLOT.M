function isoplot(filename,first_channel,second_channel,selected_channel,...
                 direction_matrix,n_directions,n_lines,sampling_rate,...
                 begin_of_signal,end_of_signal,plot_kind)

% function isoplot(filename,first_channel,second_channel,selected_channel,...
%                  direction_matrix,n_directions,n_lines,sampling_rate,...
%                  begin_of_signal,end_of_signal,plot_kind)
%
% ISOPLOT plots a 3D plot of isolevels at a selected frequency. Two different
%         kind of plots are available:
%          - Sphere plot: Combined 3D plot of the isolevels in a 3D sphere
%                         with contour lines and pseudocolors
%          - 3D-Plane plot: Combined 3D plot of the isolevels above a plane
%                           in three dimensions with contour lines and
%                           pseudocolors
%         Only one channel is used for plotting, even when the file includes
%         two channels and both channels are passed with the arguments. This
%         done for economizing memory (see MATLAB UserGuide p. 2-136).
%
%         Description of the input arguments:
%          - filename: name of the loaded file for the title of the plots
%          - first_channel: matrix with the absolute magnitude of the signals
%                           first channel of the signals. Each row represents
%                           one direction, each column one sampling line.
%          - second_channel: matrix with the absolute magnitude of the signals
%                            second channel. Each row represents one direction,
%                            each column one sampling line. If the file includes
%                            only one channel, this matrix is empty.
%          - selected_channel: channel to plot (must be '1' or '2')
%          - direction_matrix: matrix including the elevation and azimuth
%                              angle of the different directions. The first
%                              row includes the elevation, the second row
%                              the azimuth. Each column represents a
%                             direction.
%          - n_directions: number of directions
%          - n_lines: number of lines for each direction and each channel
%          - sampling_rate
%          - begin_of_signal
%          - end_of_signal
%          - plot_kind: selects the kind of 3D-plot (must be 'plane', 
%                       'sp_level' or 'sp_phase')
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


% select frequency of the isoline

phase_flag=0; % K.H
if(strcmp(plot_kind,'sp_phase')==1) phase_flag=1; end; % k.H 

isofrequency = [];
while isempty(isofrequency)
  isofrequency = input('Selected frequency: ');
  if ( (isofrequency>end_of_signal)|(isofrequency<begin_of_signal) )
    isofrequency = [];
  end;
end;

% transformation of selected frequency to corresponding line number

resolution = (end_of_signal - begin_of_signal) / (n_lines-1);
selected_line = round( (isofrequency - begin_of_signal) / resolution + 1 );
isofrequency = begin_of_signal + (selected_line - 1) * resolution;

% selected magnitudes of selected channel at the choosen line

if selected_channel==1
               magn = first_channel(:,selected_line); 
elseif selected_channel==2
           magn = second_channel(:,selected_line);
end;

% generate higher resolution of the sampled data array at the selected frequency


min_level =(min(magn));
max_level =(max(magn));


if (phase_flag==0)
  disp(['Signal mimimum level: ',num2str(min_level),' dB']);
  disp(['Signal maximum level: ',num2str(max_level),' dB']);
  dbsteps = [];
  while isempty(dbsteps)      % selecting contour line difference
    dbsteps = input('contour line difference in [dB] ? ');
     if dbsteps<0
        dbsteps = -dbsteps;
     end; 
   end;
else

  disp(['Signal mimimum phase: ',num2str(min_level),' pi']);
  disp(['Signal maximum phase: ',num2str(max_level),' pi']);
  dbsteps = [];
  while isempty(dbsteps)      % selecting contour line difference
    dbsteps = input('contour line difference in [pi] ? ');
     if dbsteps<0
        dbsteps = -dbsteps;
     end; 
   end;
end;

minline = dbsteps * fix(min_level/dbsteps);
maxline = dbsteps * fix(max_level/dbsteps);
number_of_contours = (maxline-minline)/dbsteps + 1;
level = linspace(minline,maxline,number_of_contours);

if (phase_flag==0)
[level_grid,level_color,level_contour,azi_interp,...
 elev_interp] = iso_fit(magn,direction_matrix,n_directions,min_level, max_level,level);
else
[level_grid,level_color,level_contour,azi_interp,...
 elev_interp] = iso_phfit(magn,direction_matrix,n_directions,min_level, max_level,level);
end

% plot selected kind of plot


if (strcmp(plot_kind,'plane'))          % mesh and contour plot of all levels at
                               % selected frequency
  iso_mesh(isofrequency,level_grid,level_color,level_contour,...
           azi_interp,elev_interp);
  if selected_channel==1
    title(['isolevel: first channel  ',filename,'    ',num2str(isofrequency),' Hz']);
  elseif selected_channel==2
    title(['isolevel: second channel  ',filename,'    ',num2str(isofrequency),' Hz']);
  end;

elseif (strcmp(plot_kind,'sp_level')==1)                           % plotting isolines on a sphere
  isoglobe(isofrequency,level_grid,level_color,level_contour,...
           azi_interp,elev_interp);
  if selected_channel==1
    title(['isolevel: first channel  ',filename,'    ',num2str(isofrequency),' Hz']);
     elseif selected_channel==2
    title(['isolevel: second channel  ',filename,'    ',num2str(isofrequency),' Hz']);
  end;


elseif (strcmp(plot_kind,'sp_phase')==1)
  isoglobe(isofrequency,level_grid,level_color,level_contour,...
           azi_interp,elev_interp);
  if selected_channel==1
    title(['isophase: first channel  ',filename,'    ',num2str(isofrequency),' Hz']);
  elseif selected_channel==2
    title(['isophase: second channel  ',filename,'    ',num2str(isofrequency),' Hz']);
  end;
end;

colormap(jet)
view(-90,90)
prcntrl
axis off

end;
