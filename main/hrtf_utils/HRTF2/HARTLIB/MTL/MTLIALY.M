% MTLIALY
% derive from ia_lyse

clc;
choice=0;
data_loaded = 0;
evaluated=0;
binaural=0;

while choice~=8
  choice=menu('Interaural analyse','Load data',...
              'Evaluation',...
              'Save data',...
              'Sphere-Plot of level difference',...
              'Sphere-Plot of time difference',...
              'Plot level and time difference',...
              'Level versus Time difference',...
              'Exit');
%disp(choice);

if (choice==1) 
    [filename,pathname] = uigetfile('*.mtl','Reading MTLSIG-Format');
    if (filename == 0) return; end; 
     filetype=mtltype([pathname filename]);
     if (filetype~=3)&(filetype~=11)
        disp('File has incorrect filetype !');
        return;
     end;
  
    if (filetype==3) 
     [filetype,info_blocks,...
      n_channels,n_lines, sampling_rate,...
      first_line,last_line,n_directions,...
      comment1, comment2] =mtlrh([pathname filename]);
      if ((n_channels/n_directions) ==2) binaural=1; end;
      if (binaural==0) disp('Input file has be binaural'); return; end;
   
      xaxis=linspace(first_line,last_line,n_lines/2+1);
      xaxis=xaxis(0:n_lines/2);

      % evaluate channel by channel
      for index=1:n_directions 
          [file_format, file_type, H_l,H_r, direction_matrix,...
          delay_l, delay_r, n_directions, n_channels, n_lines,...
          first_frq, last_frq]=mtlrt3([pathname filename], index);

          n_directions=1;
       
         [level_difference(index),time_difference(index),frq_axis] = ia_evalx(H_l,H_r... 
                                                                            ,n_directions,...
                                                         n_lines,sampling_rate,first_frq,...
                                                         last_frq,xaxis);

       end;
       
       data_loaded = 1;
       evaluated=1;

      end;  % filetype 3
      
      if (filetype==11)
      % read header
      % readfile
      end; % filetype 11

     % now show it
       figure();
       ia_globe(level_difference,direction_matrix,frq_axis,'level_difference');
       figure();
       ia_globe(time_difference,direction_matrix,frq_axis,'time_difference');

  




    break;
    %disp(filetype);

    end;
break;
end;

end;



