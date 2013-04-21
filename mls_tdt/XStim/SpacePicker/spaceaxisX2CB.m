function []= spaceaxisX2CB;

%spaceaxisX2CB:  button-down callback function for VAR space axis
%in spacePicker

global GUI;

spaceaxes = gcbo;
point1 = get(gca,'CurrentPoint');  %drag-select point
finalRect = rbbox;                   
point2 = get(gca,'CurrentPoint');    
point1 = point1(1,1:2);              
point2 = point2(1,1:2);   
dirmatrix = mkdirmatrix;


if GUI.cursormode == 2;
    if abs(point1(1)-point2(1))+abs(point1(2)-point2(2))<2;       %single point selection
        for i = 1:size(dirmatrix,2);           %find nearest hrtf location in dirmatrix to selected point
            distance(i) = sqrt((point1(1) - dirmatrix(2,i))^2 + (point1(2) - dirmatrix(1,i))^2);
        end
        temp = dirmatrix(:,distance==min1(distance))';
        GUI.az1 = temp(2);
        GUI.el1 = temp(1);
        if abs(GUI.az1) + abs(GUI.el1)<= 90              %selected point is within region sampled for HRTFs
            line (GUI.az1, GUI.el1, 'Color', [0.1 0.1 0.7], 'LineStyle', 'none', 'Marker','.', 'HitTest', 'off');
            if isempty(GUI.locations1)              %add point to selected locations
                GUI.locations1(1,1) = GUI.el1;
                GUI.locations1(1,2) = GUI.az1;
            else
                nLocs = size(GUI.locations1,1);
                GUI.locations1(nLocs+1, 1) = GUI.el1;
                GUI.locations1(nLocs+1, 2) = GUI.az1;
            end
        else                                     %selected point is outside of frontal 90 degrees
            line (GUI.az1, GUI.el1, 'Color', [0.75, 0.75, 1],...
                'LineStyle', 'none', 'Marker','+');
        end
    else                                                      %drag selection case
        minAz = min([point1(1), point2(1)]); 
        maxAz = max([point1(1), point2(1)]);
        minEl = min([point1(2), point2(2)]);
        maxEl = max([point1(2), point2(2)]);
        incAz = find(dirmatrix(2,:)>=minAz & dirmatrix(2,:)<=maxAz);
        incEl = find(dirmatrix(1,:)>=minEl & dirmatrix(1,:)<=maxEl);
        loc_index = intersect(incAz, incEl);
        for i = 1:length(loc_index)
            sel_locs (i, 1) = dirmatrix(1,loc_index(i));
            sel_locs (i, 2) = dirmatrix(2,loc_index(i));
        end
        if ~isempty(loc_index)
        if isempty(GUI.locations1)
            GUI.locations1 = sel_locs;
        else
            GUI.locations1 = [GUI.locations1' sel_locs']';
        end
        end
    end
    %==============================================================================
    %Delete Points: SOUND1 (cursormode = 3)
    %==============================================================================
elseif GUI.cursormode == 3
    if ~isempty(GUI.locations1)
        if abs(point1(1)-point2(1))+abs(point1(2)-point2(2))<2; %single point selection 
            dirmatrix = GUI.locations1';
            for i = 1:size(dirmatrix,2);           %find nearest hrtf location in dirmatrix to selected point
                distance(i) = sqrt((point1(1) - dirmatrix(2,i))^2 + (point1(2) - dirmatrix(1,i))^2);
            end
            near_location = find(distance==min1(distance));
            GUI.locations1(near_location,:) = [];
        else                                            % multiple point selection
            minAz = min([point1(1), point2(1)]);
            maxAz = max([point1(1), point2(1)]); 
            minEl = min([point1(2), point2(2)]);
            maxEl = max([point1(2), point2(2)]);
            incAz = find(GUI.locations1(:,2)>=minAz & GUI.locations1(:,2)<=maxAz);
            incEl = find(GUI.locations1(:,1)>=minEl & GUI.locations1(:,1)<=maxEl);
            loc_index = intersect(incAz, incEl);
            n_sel_locs = length(loc_index);  
            for i = 0:(n_sel_locs - 1);
                GUI.locations1(loc_index(n_sel_locs-i),:) = [];
            end
        end
    end
elseif GUI.cursormode == 4                % falling diagonal
        for i = 1:size(dirmatrix,2);           %find nearest hrtf location in dirmatrix to selected point
            distance(i) = sqrt((point1(1) - dirmatrix(2,i))^2 + (point1(2) - dirmatrix(1,i))^2);
        end
        temp = dirmatrix(:,distance==min1(distance))';
        for iLoc = 1:size(dirmatrix,2)
            if sum(dirmatrix(:,iLoc)) == sum(temp)
                nLocs = size(GUI.locations1,1);
                GUI.locations1(nLocs+1, 1) = dirmatrix(1,iLoc);
                GUI.locations1(nLocs+1, 2) = dirmatrix(2,iLoc);
            end
        end
elseif GUI.cursormode == 5                % rising diagonal
        for i = 1:size(dirmatrix,2);           %find nearest hrtf location in dirmatrix to selected point
            distance(i) = sqrt((point1(1) - dirmatrix(2,i))^2 + (point1(2) - dirmatrix(1,i))^2);
        end
        temp = dirmatrix(:,distance==min1(distance))';
        for iLoc = 1:size(dirmatrix,2)
            if diff(dirmatrix(:,iLoc)) == diff(temp)
                nLocs = size(GUI.locations1,1);
                GUI.locations1(nLocs+1, 1) = dirmatrix(1,iLoc);
                GUI.locations1(nLocs+1, 2) = dirmatrix(2,iLoc);
            end
        end

else disp('ERROR: unsupported cursormode')
end

SetInfo_spacePicker;