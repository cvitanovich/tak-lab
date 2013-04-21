function []= spaceaxisX2CB;

%spaceaxisX2CB:  button-down callback function for VAR space axis
%in spaceConsole

global GUI;

spaceaxes = gcbo;
point1 = get(gca,'CurrentPoint');  %drag-select point
finalRect = rbbox;                   
point2 = get(gca,'CurrentPoint');    
point1 = point1(1,1:2);              
point2 = point2(1,1:2);   

if GUI.cursormode == 2;
    az_sel = point1(1);
    el_sel = point1(2);
    dirmatrix = mkdirmatrix(GUI.spaceres);
    if abs(point1(1)-point2(1))+abs(point1(2)-point2(2))<2; %single point selection
        nlocs = size(dirmatrix,2);
        distance = zeros(1, nlocs);
        for i = 1:nlocs; %find nearest hrtf location in dirmatrix to selected point
            distance(i) = sqrt((az_sel - dirmatrix(2,i))^2 + (el_sel - dirmatrix(1,i))^2);
        end
        mindist = min(distance);
        near_location = find(distance==mindist);
        if length(near_location)>1
            near_location = near_location(1);
        end
        GUI.az1 = dirmatrix(2, near_location);
        GUI.el1 = dirmatrix(1, near_location);
        if abs(GUI.az1) + abs(GUI.el1)<= 90              %selected point is within region sampled for HRTFs
            temp = findobj(gcf,'Tag','azimuth');
            set(temp,'string', num2str(GUI.az1));
            temp = findobj(gcf,'Tag','elevation');
            set(temp,'string', num2str(GUI.el1));
            line (GUI.az1, GUI.el1, 'Color', [0.1 0.1 0.7], 'LineStyle', 'none', 'Marker','.', 'HitTest', 'off');
            if isempty(GUI.locations1)==1;              %add point to selected locations
                GUI.locations1(1,1) = GUI.el1;
                GUI.locations1(2,1) = GUI.az1;
            else
                ind = find(GUI.locations1(1,:)==GUI.el1 & GUI.locations1(2,:)==GUI.az1);
                if isempty(ind)
                    nlocations = size(GUI.locations1,2);
                    GUI.locations1(1,nlocations+1) = GUI.el1;
                    GUI.locations1(2,nlocations+1) = GUI.az1;
                end
            end
        else                                     %selected point is outside of frontal 90 degrees
            line (GUI.az1, GUI.el1, 'Color', [0.75, 0.75, 1],...
                'LineStyle', 'none', 'Marker','+');
        end
    else   %drag selection case
        minAz = min([point1(1), point2(1)]); 
        maxAz = max([point1(1), point2(1)]);
        minEl = min([point1(2), point2(2)]);
        maxEl = max([point1(2), point2(2)]);
        incAz = find(dirmatrix(2,:)>=minAz & dirmatrix(2,:)<=maxAz);
        incEl = find(dirmatrix(1,:)>=minEl & dirmatrix(1,:)<=maxEl);
        loc_index = intersect(incAz, incEl);
        sel_locs = zeros(2,length(loc_index));
        for i = 1:length(loc_index);
            sel_locs (1,i) = dirmatrix(1,loc_index(i));
            sel_locs (2,i) = dirmatrix(2,loc_index(i));
        end
        if isempty(GUI.locations1);
            GUI.locations1 = sel_locs;
        else
            for i = 1:length(loc_index)
                ind = find(GUI.locations1(1,:)==sel_locs(1,i) & GUI.locations1(2,:)==sel_locs(2,i));
                if isempty(ind)
                    GUI.locations1 = [GUI.locations1 sel_locs(:,i)];
                end
            end
        end
        cla;
        refresh_spaceaxis;
    end
    %==============================================================================
    %Delete Points: SOUND1 (cursormode = 3)
    %==============================================================================
elseif GUI.cursormode == 3;
    if ~isempty(GUI.locations1)
        if abs(point1(1)-point2(1))+abs(point1(2)-point2(2))<2; %single point selection                                                         
            az_sel = point1(1);
            el_sel = point1(2);
            nlocs = size(GUI.locations1,2);
            distance = zeros(1, nlocs);
            for i = 1:nlocs;
                distance(i) = sqrt((az_sel - GUI.locations1(2,i))^2 + (el_sel - GUI.locations1(1,i))^2);
            end
            mindist = min(distance);
            near_location = find(distance==mindist);
            if length(near_location)>1
                near_location = near_location(1);
            end
            GUI.locations1(:,near_location) = [];
            refresh_spaceaxis;
        else   %multiple point selection
            minAz = min([point1(1), point2(1)]);
            maxAz = max([point1(1), point2(1)]); 
            minEl = min([point1(2), point2(2)]);
            maxEl = max([point1(2), point2(2)]);
            incAz = find(GUI.locations1(2,:)>=minAz & GUI.locations1(2,:)<=maxAz);
            incEl = find(GUI.locations1(1,:)>=minEl & GUI.locations1(1,:)<=maxEl);
            loc_index = intersect(incAz, incEl);
            n_sel_locs = length(loc_index);  
            for i = 0:(n_sel_locs - 1);
                GUI.locations1(:,loc_index(n_sel_locs-i)) = [];
            end
            refresh_spaceaxis;
        end
    end
else disp('ERROR: unsupported cursormode')
end

SetInfo_spacePicker;