function [] = SetInfo_spacePicker;

% SetInfo_spacePicker
global GUI
global H
global FN
global REC_INFO

if exist1('H.pickerfig') figure(H.pickerfig);	end

% check for UseLastLocations_flag
if get(H.UseLastLocations,'value')
    GUI.UseLastLocations_flag =1;
    eval(['load ' FN.current_path 'Locations_current;'])
    GUI.locations1 = locations;
    clear locations
else
    GUI.UseLastLocations_flag =0;
    %save locations to file
    if ~isempty(GUI.locations1)
        locations = GUI.locations1;
        eval(['save ' FN.current_path 'Locations_current locations;'])
        clear locations
    end
end

% all locations PB
if get(H.allPB,'value')
    GUI.locations1 = mkdirmatrix(GUI.spaceres);
end 
 
% addpoint
if get(H.addpoint,'value')
    GUI.cursormode = 2;
    set(findobj(gcf,'Tag','mouse_point1'),'value',0);
    set(findobj(gcf,'Tag','delete_point1'),'value',0);
end

% delete point
if get(H.delpoint,'value')
    GUI.cursormode = 3;
    set(findobj(gcf,'Tag','mouse_point1'),'value', 0);
    set(findobj(gcf,'Tag','add_point1'),'value', 0);
end

% resolution
switch get(H.resolution, 'value')
    case 1
        GUI.spaceres = 'max';
    case 2
        GUI.spaceres = '5deg';
    case 3
        GUI.spaceres = '10degOffset';
    case 4 
        GUI.spaceres = '10deg';
    case 5
        GUI.spaceres = '20degOffset';
    case 6
        GUI.spaceres = '20deg';
    case 7
        GUI.spaceres = '9Locs';
        iloc = 0;
        for el = 30:-30:-30
            for az = -30:30:30
                iloc = iloc+1;
                GUI.locations1(:,iloc) = [el;az];
            end
        end
end

%% clear PB
if get(H.clearPB,'value')
    GUI.locations1 = [];
end

% remove edges
if get(H.removeEdges,'value') & ~isempty(GUI.locations1)
    ind = find(abs(GUI.locations1(1,:)) + abs(GUI.locations1(2,:)) < 90);
    GUI.locations1 = GUI.locations1(:,ind);
end

% update display
refresh_spaceaxis
