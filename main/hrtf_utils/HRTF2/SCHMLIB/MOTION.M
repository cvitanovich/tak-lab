% MOTION shows a motion consisting of all directions, included in the matrix
%        'channel'. 'channel' must be provided outside this M-File, and
%        includes the signals to plot in a motion. Depending on the status-
%        variable 'motion_dim' a 2D or 3D motion is produced. The motion is
%        layouted with a title, axes ranges and axes labels. A pushbutton
%        'Replay' allows the repetition of the motion.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

if filetype~=2
  channel = abs2db(channel);
end;

clf;
axis('off');
m = moviein(n_directions);
if motion_dim==2    % 2D - motion
 axis('on');
 if scaling=='lin'
   for j=1:n_directions
     signal = channel(j,:);
     plot(xaxis,signal);
     if filetype==2           % time signal
      if channel_number==1
       t_layout(file_format,filetype,filename,1,'impulse response first channel',j,...
                direction_matrix,xaxis,n_lines);
      elseif channel_number==2
       t_layout(file_format,filetype,filename,1,'impulse response second channel',j,...
                direction_matrix,xaxis,n_lines);
      end;
     else                     % frequency signal
       if channel_number==1
         m_layout(file_format,filetype,filename,1,'transfer function first channel',j,...
                  direction_matrix,xaxis,n_lines);
       elseif channel_number==2
         m_layout(file_format,filetype,filename,1,'transfer function second channel',j,...
                  direction_matrix,xaxis,n_lines);
       end;
     end;
     m(:,j) = getframe;
   end;
 elseif scaling=='log'
   for j=1:n_directions
     signal = channel(j,:);
     semilogx(xaxis,signal);
     if filetype==2           % time signal
      if channel_number==1
       t_layout(file_format,filetype,filename,1,'impulse response first channel',j,...
                direction_matrix,xaxis,n_lines);
      elseif channel_number==2
       t_layout(file_format,filetype,filename,1,'impulse response second channel',j,...
                direction_matrix,xaxis,n_lines);
      end;
     else                     % frequency signal
       if channel_number==1
         m_layout(file_format,filetype,filename,1,'transfer function first channel',j,...
                  direction_matrix,xaxis,n_lines);
       elseif channel_number==2
         m_layout(file_format,filetype,filename,1,'transfer function second channel',j,...
                  direction_matrix,xaxis,n_lines);
       end;
     end;
     m(:,j) = getframe;
   end;
 end;
else                % 3D - motion
 axis('on');
 for j=1:n_directions
   basic = ones(size(channel));
   signal = channel(j,:);
   basic(j,:) = signal;
  
   [yaxis,lab,ang] = y_axis(direction_matrix,n_directions);
   mesh(xaxis,yaxis,basic);
   grid on;
   z_status = 0;
   mesh_lay;
   if filetype==2           % time signal
    if channel_number==1
     title(['impulse response first channel: ',filename]);
    elseif channel_number==2
     title(['impulse response second channel: ',filename]);
    end;
   else                     % frequency signal
    if channel_number==1
     title(['transfer function first channel: ',filename]);
    elseif channel_number==2
     title(['transfer function second channel: ',filename]);
    end;
   end;
   m(:,j) = getframe;
 end;
end;

movie(m,1,fps);

% Push button

callback = [ ...
   'set(uip,''back'',[1 2/3 1/3]);' ...
   'movie(m,1,fps); ' ...
   'set(uip,''back'',''default'');' ];

uip = uicontrol('style','push','pos',[20 10 60 20],'string','Replay',...
   'call',callback);
