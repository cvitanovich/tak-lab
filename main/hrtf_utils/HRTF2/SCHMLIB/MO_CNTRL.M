% MO_CNTRL is a submenu called from 'analyse', for the control of motions.
%          The motion displays always the signals magnitude. Depending on
%          the kind of motion and the number of channels, different menus
%          are provided. Motions in 2D and 3D are possible.
%
%          The following choices are available in the different menus:
%           - Speed: changes the motion speed; default is 3 frames per
%                    second
%           - Change to 2D / 3D: switches between the different menus for 
%                                2D or 3D motions
%           - Close: closes motion menu and returns to the menu provided
%                    from 'analyse'
%           - Play: start of the motion when the signals have only one
%                   channel
%           - Channel 1 / 2: motion start of the first or second
%                            channel when the signals have two channels
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

fps = 3;    %frames per second
motion_dim = 2;         % Standard : 2D-motion
figure;

if n_channels==1

  subchoice=0;
  while subchoice~=4
     if motion_dim==2

       subchoice=menu('2D - Movie Control',...
                      'Speed',...
                      'Play',...
                      'Change to 3D',...
                      'Close');
     else
       subchoice=menu('3D - Movie Control',...
                      'Speed',...
                      'Play',...
                      'Change to 2D',...
                      'Close');
     end;


     if subchoice==1
       movie_seq = mo_speed;
     elseif subchoice==2
       channel = first_channel;
       channel_number = 1;
       motion;
     elseif subchoice==3
        if motion_dim==2
          motion_dim=3;
        else
          motion_dim=2;
        end;       
     elseif subchoice==4
        delete(gcf);
        break;
     end;
  end;

elseif n_channels==2

  subchoice=0;
  while subchoice~=5
     if motion_dim==2
      subchoice=menu('2D - Movie Control',...
                     'Speed',...
                     'Channel 1',...
                     'Channel 2',...
                     'Change to 3D',...
                     'Close');
     else
      subchoice=menu('3D - Movie Control',...
                     'Speed',...
                     'Channel 1',...
                     'Channel 2',...
                     'Change to 2D',...
                     'Close');
     end;  
    
     if subchoice==1
       fps = mo_speed;
     elseif subchoice==2
       channel = first_channel;
       channel_number = 1;
       motion;
     elseif subchoice==3
       channel = second_channel;
       channel_number = 2;
       motion;
     elseif subchoice==4
       if motion_dim==2
          motion_dim=3;
       else
          motion_dim=2;
       end;      
     elseif subchoice==5
        delete(gcf);
        break;
     end;
  end;

end;