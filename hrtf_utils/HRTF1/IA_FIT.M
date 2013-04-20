function [level_grid,level_color,level_contour,azi_interp,...
          elev_interp] = ia_fit(difference,direction_matrix,n_directions,type)

% function [level_grid,level_color,level_contour,azi_interp,...
%           elev_interp] = ia_fit(difference,direction_matrix,n_directions,type)
%
% IA_FIT linear interpolation of the interaural level or time differences of
%        the measured directions to directions, that are not measured, to
%        get a better resolution.
%        All measured directions of a certain elevation angle are linear
%        interpolated about an range of 360— in steps of 3—. This is done
%        for all different elevation angles between the minimum and the
%        maximum elevation. The interpolated azimuth angles are summarized
%        to a matrix that, then is interpolated about the elevation angles
%        with steps of 2—.
%        The resolution of the measured directions must be an integer multiple
%        of 3 for the azimuth and an integer multiple of 2 for the elevation,
%        or the measured directions are not included in the output argument
%        'level_grid'. 
%        The interpolated interaural differences are used for detecting
%        contour lines and assigning colors to the levels.         
%
%        Description of the input arguments:
%        - difference: vector including the interaural time or level differences
%                      of the measured directions at a selected frequency
%        - direction_matrix: matrix including the elevation and azimuth
%                            angle of the different directions. The first
%                            row includes the elevation, the second row
%                            the azimuth. Each column represents a
%                            direction.
%        - n_directions: number of directions
%        - type: separates between interaural 'level_difference' and 'time_difference'
%                for interpolation
%
%         Description of the output arguments:
%          - level_grid: matrix including the interpolated interaural differences
%                        Each column represents a certain azimuth angle,
%                        each row represents a certain elevation angle.
%          - level_color: matrix including the pseudocolors belonging to the
%                         interpolated interaural differences. Each column
%                         represents a certain azimuth angle, each row
%                         represents a certain elevation angle.
%          - level_contour: matrix including the output arguments of the function
%                           'contourc'
%          - azi_interp: vector including the azimuth angles
%          - elev_interp: vector including the elevation angles
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter


% set number of contour lines

if strcmp(type,'level_difference')
  min_diff = min(difference);
  max_diff = max(difference);
  disp(['Mimimum difference: ',num2str(min_diff), 'dB']);
  disp(['Maximum difference: ',num2str(max_diff),' dB']);
elseif strcmp(type,'time_difference')
  difference = difference * 10^6;
  min_diff = min(difference);
  max_diff = max(difference);
  disp(['Mimimum difference: ',num2str(min_diff),' * 10^(-6) sec']);
  disp(['Maximum difference: ',num2str(max_diff),' * 10^(-6) sec']);
end;

diffsteps = [];
while isempty(diffsteps)      % selecting contour line difference
  if strcmp(type,'level_difference')
    diffsteps = input('contour line difference in [dB] ? ');
  elseif strcmp(type,'time_difference')
    diffsteps = input('contour line difference in [10^(-6) sec] ? ');
  end;
  if diffsteps<0
    diffsteps = -diffsteps;
  end; 
end;
minline = diffsteps * fix(min_diff/diffsteps);
maxline = diffsteps * fix(max_diff/diffsteps);
number_of_contours = (maxline-minline)/diffsteps + 1;
level = linspace(minline,maxline,number_of_contours);

% set spatial borders for approximation

min_elev = min(direction_matrix(1,:));
max_elev = max(direction_matrix(1,:));
min_azi = min(direction_matrix(2,:));
max_azi = max(direction_matrix(2,:));

delta=abs(max_azi-min_azi);
if (delta>180) delta=abs(delta-360); end;% added by K. Hartung 19.2.95
%if (delta<=45) min_azi=-max_azi; end; % added k. Hartung 19.2.95
if (delta<=45) max_azi=min_azi+360; end;% added k. Hartung 19.2.95


% detect number of different elevation angles

number_of_elevations = 0;
rest_directions = sort(direction_matrix(1,:));
while isempty(rest_directions)==0
  number_of_elevations = number_of_elevations + 1;
  elev_array(number_of_elevations) = rest_directions(1);
  rest_directions = rest_directions(find(rest_directions~=rest_directions(1)));
end;

% detect number of different azimuth angles

number_of_azimuths = 0;
rest_directions = sort(direction_matrix(2,:));
while isempty(rest_directions)==0
  number_of_azimuths = number_of_azimuths + 1;
  azi_array(number_of_azimuths) = rest_directions(1);
  rest_directions = rest_directions(find(rest_directions~=rest_directions(1)));
end;

% Interpolate the magnitude (abs) of a certain elevation angle for all azimuth angles.
% The first and last azimuth angle are added with a angle shift of +-360 degrees
% for linear interpolation. Then resulting vectors are inserted in a matrix, with 
% the elevation angles as rows and the azimuth angles as columns. This matrix is
% interpolated now about the elevation angles.

magn_interp = [];

disp('Interpolation of Isopoints');
for i=1:number_of_elevations
  elev = elev_array(i);
  azi_index = find(elev==direction_matrix(1,:));
  azi_amount = length(azi_index);
  orig_azi = direction_matrix(2,azi_index);
 

 [orig_azi,sort_index] = sort(orig_azi);  % added by K. Hartung 13.2.95
  azi_index=azi_index(sort_index);        % added by K. Hartung 13.2.95


  
  orig_diff = difference(azi_index);
  approx_diff = [azi_amount+2];
  approx_azi = [azi_amount+2];
  approx_azi(1) = orig_azi(azi_amount)-360;
  approx_diff(1) = orig_diff(azi_amount);
  approx_azi(2:azi_amount+1) = orig_azi;
  approx_diff(2:azi_amount+1) = orig_diff;
  approx_azi(azi_amount+2) = orig_azi(1)+360;
  approx_diff(azi_amount+2) = orig_diff(1);
%pause
%  azi_interp = -180:3:180;    % interpolate in 3— steps
 azi_interp = min_azi:3:max_azi;    % interpolate in 3— steps

 level_interp(i,:) = (interp1(approx_azi,approx_diff,azi_interp))';
end;

elev_interp = min_elev:2:max_elev;     % interpolate in 2— steps
level_grid = interp2(azi_interp,elev_array,level_interp,azi_interp,elev_interp);

% evaluate color levels for the level_grid

disp('Evaluate surface colors');
fig_h = figure;
surf_h = surface(azi_interp,elev_interp,level_grid);
level_color = get(surf_h,'Cdata');
set(surf_h,'visible','off');

close(fig_h);

% set contour lines

disp('Evaluate contour lines');
level_contour = contourc(azi_interp,elev_interp,level_grid,level);

