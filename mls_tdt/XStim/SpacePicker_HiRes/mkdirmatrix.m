function [dir_matrix] = mkdirmatrix(spaceres);

%[dir_matrix] = mkdirmatrix(spaceres); generates matrix of spatial locations in double polar coordinates
%dir_matrix: column1 = elevations, column2 = azimuths
%argin spaceres: determines desired spatial resolution

global GUI
global XStimParams
global FN

if strcmp('max',GUI.spaceres);
    inc = 2;
else
    inc = 5;
end

iloc = 0;
for el = -90:inc:90
    for az = -90:inc:90
        if (abs(el)+abs(az)) <= 90
            iloc = iloc+1;
            temp_dir(:,iloc) = [el; az];
        end
    end
end
Ndirs = iloc;    
    
if strcmp('max',GUI.spaceres);
    dir_matrix =  temp_dir;
    
 elseif strcmp('5deg',GUI.spaceres)
    n = 0;
    for i = 1:Ndirs
        if (rem(temp_dir(1,i),5)==0) & (rem(temp_dir(2,i),5)==0)
            n = n+1;
            dir_matrix(:,n) = temp_dir(:,i);
        end
    end
    
elseif strcmp('10degOffset',GUI.spaceres);
    n = 0;
    for i = 1:Ndirs
        if (rem((temp_dir(1,i) + temp_dir(2,i)),10) == 0) & (rem(temp_dir(1,i),5)==0) & (rem(temp_dir(2,i),5)==0)
            n = n+1;
            dir_matrix(:,n) = temp_dir(:,i);
        end
    end
elseif strcmp('10deg',GUI.spaceres)
    n = 0;
    for i = 1:Ndirs
        if (rem(temp_dir(1,i),10)==0) & (rem(temp_dir(2,i),10)==0)
            n = n+1;
            dir_matrix(:,n) = temp_dir(:,i);
        end
    end
    
elseif strcmp('20degOffset',GUI.spaceres)
    n = 0;
    for i = 1:Ndirs
        if (rem((temp_dir(1,i) + temp_dir(2,i)-10),20) == 0) & (rem(temp_dir(1,i),10)==0) & (rem(temp_dir(2,i),10)==0)
            n = n+1;
            dir_matrix(:,n) = temp_dir(:,i);
        end
    end

    
elseif strcmp('20deg',GUI.spaceres)
    n = 0;
    for i = 1:Ndirs
        if (rem(temp_dir(1,i),20)==0) & (rem(temp_dir(2,i)-10,20)==0)
            n = n+1;
            dir_matrix(:,n) = temp_dir(:,i);
        end
    end

elseif strcmp('9Locs',GUI.spaceres)
    iloc = 0;
    for el = 30:-30:-30
        for az = 30:-30:30
            iloc = iloc+1;
            dir_matrix(:,iloc) = [el;az];
        end
    end
    
else disp ('ERROR: requested unsupported spatial resolution')
end


