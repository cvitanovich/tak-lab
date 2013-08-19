function [] = SetInfo_spacePicker_ted;

% SetInfo_spacePicker
global GUI
global H
global FN
global XStimParams

% check for UseLastLocations_flag
if get(H.UseLastLocations, 'value')
    if H.shape_flag == 1;
		GUI.UseLastLocations_flag = 1;
		eval(['load ' FN.current_path 'Locations_current1;'])
		GUI.locations1 = locations;
		for count1 = 1:size(GUI.locations1, 1)
			for count2 = (count1 + 1):size(GUI.locations1, 1)
				if GUI.locations1(count1, 1) == GUI.locations1(count2, 1) & GUI.locations1(count1, 2) == GUI.locations1(count2, 2)
					GUI.locations1(count1, :) = [];
					break;
				end
			end
		end
		clear locations
	else
		GUI.UseLastLocations_flag = 1;
		eval(['load ' FN.current_path 'Locations_current2;'])
		GUI.locations1 = locations;
		for count1 = 1:size(GUI.locations1, 1)
			for count2 = (count1 + 1):size(GUI.locations1, 1)
				if GUI.locations1(count1, 1) == GUI.locations1(count2, 1) & GUI.locations1(count1, 2) == GUI.locations1(count2, 2)
					GUI.locations1(count1, :) = [];
					break;
				end
			end
		end
		clear locations
	end
else
    %save locations to file
    if ~isempty(GUI.locations1)
		locations = GUI.locations1;
		for count1 = 1:size(locations, 1)
			for count2 = (count1 + 1):size(locations, 1)
				if locations(count1, 1) == locations(count2, 1) & locations(count1, 2) == locations(count2, 2)
					locations(count1, :) = [];
					break;
				end
			end
		end
		if H.shape_flag == 1;
			eval(['save ' FN.current_path 'Locations_current1 locations;'])
			clear locations
		else
			eval(['save ' FN.current_path 'Locations_current2 locations;'])
			clear locations
		end
	end
end

GUI.UseLastLocations_flag = 0;

% remove edges
if get(H.removeEdges,'value') & ~isempty(GUI.locations1)
    ind = find(abs(GUI.locations1(:,1)) + abs(GUI.locations1(:,2)) < 90);
    GUI.locations1 = GUI.locations1(ind,:);
end

% set cursormode for add or delete pts
if get(H.add_point1,'value')
    GUI.cursormode = 2;
end
if get(H.delete_point1,'Value')
    GUI.cursormode = 3;
end
if get(H.falling_diagonal,'Value')
    GUI.cursormode = 4;
elseif get(H.rising_diagonal,'Value')
    GUI.cursormode = 5;
end
switch GUI.cursormode
    case 2
    set(H.delete_point1,'value', 0);
    set(H.add_point1,'BackgroundColor', 'green')
    set(H.delete_point1,'BackgroundColor',[.7 .7 .7]);
    set(H.rising_diagonal, 'BackGroundColor',[.9 .9 .9]);
    set(H.falling_diagonal, 'BackGroundColor',[.9 .9 .9]);
    case 3
    set(H.add_point1,'value', 0);
    set(H.delete_point1,'BackgroundColor', 'green')
    set(H.add_point1,'BackgroundColor', [.7 .7 .7]);
    set(H.rising_diagonal, 'BackGroundColor',[.9 .9 .9]);
    set(H.falling_diagonal, 'BackGroundColor',[.9 .9 .9]);
    case 4
    set(H.falling_diagonal, 'BackGroundColor','green')
    set(H.rising_diagonal, 'BackGroundColor',[.9 .9 .9]);
    set(H.add_point1,'BackgroundColor', 'green')
    set(H.delete_point1,'BackgroundColor',[.7 .7 .7]);
    case 5
    set(H.rising_diagonal, 'BackGroundColor','green')
    set(H.falling_diagonal, 'BackGroundColor',[.9 .9 .9]);
    set(H.add_point1,'BackgroundColor', 'green')
    set(H.delete_point1,'BackgroundColor',[.7 .7 .7]);
    otherwise
        disp('ERROR in spacepicker: incorrect cursormode');   return;
end
% select all locations for which there are hrtfs
if get(H.AllSpace,'Value')
    dirmatrix = mkdirmatrix;
    nLocs = size(dirmatrix,2);
    GUI.locations1 = zeros(nLocs, 2);
    for i = 1:nLocs;
        GUI.locations1(i,1) = dirmatrix(1,i);
        GUI.locations1(i,2) = dirmatrix(2,i);
    end
end

% clear all points
if get(H.clear,'Value')
    if GUI.cursormode < 1
        disp('ERROR: UNSUPPORTED CURSORMODE')
    elseif GUI.cursormode <= 3;
        GUI.locations1 = [];
    elseif GUI.cursormode <=6;
        GUI.locations1 = [];
    else
        disp('ERROR: UNSUPPORTED CURSORMODE')
    end
end

% check resolution
str = get(H.resolution,'String');
val = get(H.resolution,'Value');
if ~strcmp(XStimParams.GUIspaceres, deblank(str(val,:)))
switch val
    case 1
        XStimParams.GUIspaceres = 'max';
    case 2
        XStimParams.GUIspaceres = '10degOffset';
    case 3
        XStimParams.GUIspaceres = '10deg';
    case 4
        XStimParams.GUIspaceres = '20degOffset';
    case 5
        XStimParams.GUIspaceres = '20deg';
end
end

% make dirmatrix from which to choose locations
mkdirmatrix;

% update display
refresh_spaceaxis

