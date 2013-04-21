function [] = SetInfo_spacePicker;

% SetInfo_spacePicker
global GUI
global H
global FN
global XStimParams

% check for UseLastLocations_flag
if get(H.UseLastLocations,'value')
    GUI.UseLastLocations_flag =1;
    eval(['load ' FN.current_path 'Locations_current;'])
    GUI.locations1 = locations;
    clear locations
else
    %save locations to file
    if ~isempty(GUI.locations1)
        locations = GUI.locations1;
        eval(['save ' FN.current_path 'Locations_current locations;'])
        clear locations
    end
end
GUI.UseLastLocations_flag =0;

% remove edges
if get(H.picker_removeEdges,'value') & ~isempty(GUI.locations1)
    ind = find(abs(GUI.locations1(:,1)) + abs(GUI.locations1(:,2)) < 90);
    GUI.locations1 = GUI.locations1(ind,:);
end

% set cursormode for add or delete pts
if get(H.picker_add_point1,'value')
    GUI.cursormode = 2;
end
if get(H.picker_delete_point1,'Value')
    GUI.cursormode = 3;
end
if get(H.picker_falling_diagonal,'Value')
    GUI.cursormode = 4;
elseif get(H.picker_rising_diagonal,'Value')
    GUI.cursormode = 5;
end
switch GUI.cursormode
    case 2
        set(H.picker_delete_point1,'value', 0);
        set(H.picker_add_point1,'BackgroundColor', 'green')
        set(H.picker_delete_point1,'BackgroundColor',[.7 .7 .7]);
        set(H.picker_rising_diagonal, 'BackGroundColor',[.9 .9 .9]);
        set(H.picker_falling_diagonal, 'BackGroundColor',[.9 .9 .9]);
    case 3
        set(H.picker_add_point1,'value', 0);
        set(H.picker_delete_point1,'BackgroundColor', 'green')
        set(H.picker_add_point1,'BackgroundColor', [.7 .7 .7]);
        set(H.picker_rising_diagonal, 'BackGroundColor',[.9 .9 .9]);
        set(H.picker_falling_diagonal, 'BackGroundColor',[.9 .9 .9]);
    case 4
        set(H.picker_falling_diagonal, 'BackGroundColor','green')
        set(H.picker_rising_diagonal, 'BackGroundColor',[.9 .9 .9]);
        set(H.picker_add_point1,'BackgroundColor', 'green')
        set(H.picker_delete_point1,'BackgroundColor',[.7 .7 .7]);
    case 5
        set(H.picker_rising_diagonal, 'BackGroundColor','green')
        set(H.picker_falling_diagonal, 'BackGroundColor',[.9 .9 .9]);
        set(H.picker_add_point1,'BackgroundColor', 'green')
        set(H.picker_delete_point1,'BackgroundColor',[.7 .7 .7]);
    otherwise
        disp('ERROR in spacepicker: incorrect cursormode');   return;
end
% select all locations for which there are hrtfs
if get(H.picker_AllSpace,'Value')
    dirmatrix = mkdirmatrix;
    nLocs = size(dirmatrix,2);
    GUI.locations1 = zeros(nLocs, 2);
    for i = 1:nLocs;
        GUI.locations1(i,1) = dirmatrix(1,i);
        GUI.locations1(i,2) = dirmatrix(2,i);
    end
end

% clear all points
if get(H.picker_clear,'Value')
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
str = get(H.picker_resolution,'String');
val = get(H.picker_resolution,'Value');
if ~strcmp(XStimParams.GUIspaceres, deblank(str(val,:)))
    switch val
        case 1
            XStimParams.GUIspaceres = 'max';
        case 2
            XStimParams.GUIspaceres = '5deg';
        case 3
            XStimParams.GUIspaceres = '10degOffset';
        case 4
            XStimParams.GUIspaceres = '10deg';
        case 5
            XStimParams.GUIspaceres = '20degOffset';
        case 6
            XStimParams.GUIspaceres = '20deg';
    end
end

% make dirmatrix from which to choose locations
mkdirmatrix;

% update display
refresh_spaceaxis