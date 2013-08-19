function []= refresh_spaceaxis

%function []= refresh_spaceaxis
% clears space axis in SpacePicker & updates plotting of selected locations.

global GUI

%make spaceaxis current axis
if isempty(findobj('Tag','spaceaxes'))
    return
else
    axes(findobj('Tag','spaceaxes'));
end

%clear and draw patches
cla;
Xl = [-90 -90 0];
Yup = [0 90 90];
Ydn = [0 -90 -90];
Xr = [90 90 0]; 
patch(Xl, Yup, [0.9 0.9 0.9])
patch(Xl, Ydn, [0.9 0.9 0.9])
patch(Xr, Yup, [0.9 0.9 0.9])
patch(Xr, Ydn, [0.9 0.9 0.9])
if GUI.recordmode == 1;
    line([-90 0],[0 90],'Color', 'r', 'LineWidth', 1);
    line([0 90],[90 0],'Color', 'r', 'LineWidth', 1);
    line([0 90],[-90 0],'Color', 'r', 'LineWidth', 1);
    line([-90 0],[0 -90],'Color', 'r', 'LineWidth', 1);
end

% check for multi-picked locations
nLocs = size(GUI.locations1,1);
iLoc = 1;
while iLoc <= nLocs
    if ismember(GUI.locations1(iLoc,:),GUI.locations1((iLoc+1):nLocs,:),'rows');
        GUI.locations1 = GUI.locations1(setdiff(1:nLocs,iLoc),:);
        nLocs = nLocs-1;
    else
        iLoc = iLoc+1;
    end
end

%plot selected locations
if ~isempty(GUI.locations1)
    line (GUI.locations1(:,2), GUI.locations1(:,1), 'Color', [0.1 0.1 0.7], 'LineStyle', 'none',...
        'Marker','.','HitTest', 'off');
end
if ~isempty(GUI.locations2)
    line (GUI.locations2(:,2), GUI.locations2(:,1), 'Color', [0.7 0.1 0.1], 'LineStyle', 'none',...
        'Marker','.','HitTest', 'off');
end

