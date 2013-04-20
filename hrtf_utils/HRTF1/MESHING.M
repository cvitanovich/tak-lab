% MESHING Meshplot of the actual signal matrix. Sliders for changing the
%         point of view are provided, too.
% 
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


resample_factor = 10;
fig_h = figure;

disp(['IN MESHING.M JUST SENDING OUT TO y axis']);
disp(direction_matrix);
disp(n_directions);

if file_format=='MTLSIG'
  [yaxis,lab,ang] = y_axis(direction_matrix,n_directions);
elseif file_format=='SIGMAP'
  yaxis = linspace(1,n_directions,n_directions);
end;

[u_s_signal,new_xaxis,newlines] = u_sample(signal,xaxis,n_directions,...
                                           n_lines,resample_factor);
mesh(new_xaxis,yaxis,u_s_signal);
if scaling=='log'
  set(gca,'xscale','log');
end;

grid on;

mesh_lay;  % layout of the figure

v3dcntrl(fig_h);
