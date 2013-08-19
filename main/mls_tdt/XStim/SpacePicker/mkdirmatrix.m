function [dir_matrix] = mkdirmatrix();

%[dir_matrix] = mkdirmatrix(); generates matrix of spatial locations in double polar coordinates
%dir_matrix: column1 = elevations, column2 = azimuths

global XStimParams
global FN

% added Dec 12, 2005
switch XStimParams.test_type
    case 'AltIR'
        HRTF_type = 2;      %'std'
        suffix = '_std';
    case 'Space'
        HRTF_type = 2;      %'std'
        suffix = '_std';
    case 'Space2'
        HRTF_type = 2;      %'std'
        suffix = '_std';
    case '2Source'
        HRTF_type = 2;      %'std'
        suffix = '_std';
    case 'Space3'
        HRTF_type = 1;      %'eq'
        suffix = '_eq';
    case 'Mc2Source' 
        HRTF_type = 2;      %'eq'
        suffix = '_std';
    case 'Mc2SourceReversed' 
        HRTF_type = 2;      %'eq'
        suffix = '_std';

    otherwise
        HRTF_type = 1;      % default
        suffix = '_eq';
end


dir = 0;
if XStimParams.ildalone_flag 
    if FN.HRTFfiletype(2,HRTF_type) == 1
        eval(['temp_dir =  sph2dbl(mtlrdir( ''' [FN.ILA_path eval(['FN.ildalone' suffix])] '''));'])
    elseif FN.HRTFfiletype(2,HRTF_type) == 2
        eval(['load -mat ' FN.ILA_path eval(['FN.ildalone' suffix])]);
        temp_dir = dir;
    end
end
if XStimParams.itdalone_flag
    if FN.HRTFfiletype(3,HRTF_type) == 1
        eval(['temp_dir =  sph2dbl(mtlrdir( ''' [FN.ITA_path eval(['FN.itdalone' suffix])] '''));'])
    elseif FN.HRTFfiletype(3,HRTF_type) == 2
        eval(['load -mat ' FN.ITA_path eval(['FN.itdalone' suffix])]);
        temp_dir = dir;
    end
end
if XStimParams.space_flag | XStimParams.ABLalone_flag
    if FN.HRTFfiletype(1,HRTF_type) == 1
        eval(['temp_dir =  sph2dbl(mtlrdir( ''' [FN.space_path eval(['FN.space' suffix])] '''));'])
    elseif FN.HRTFfiletype(1,HRTF_type) == 2
        eval(['load -mat ' FN.space_path eval(['FN.space' suffix])]);
        temp_dir = dir;
    end
end
if XStimParams.ABLequal_flag
    if FN.HRTFfiletype(7,HRTF_type) == 1
        eval(['temp_dir =  sph2dbl(mtlrdir( ''' [FN.space_path eval(['FN.ablequal' suffix])] '''));'])
    elseif FN.HRTFfiletype(7,HRTF_type) == 2
        eval(['load -mat ' FN.space_path eval(['FN.ablequal' suffix])]);
        temp_dir = dir;
    end
end

Ndirs = size(temp_dir,2);

n = 0;
switch XStimParams.GUIspaceres
    case 'max'
        dir_matrix =  temp_dir;
        
    case '5deg'
        for i = 1:Ndirs
            if (rem(temp_dir(1,i),5)==0) & (rem(temp_dir(2,i),5)==0)
                n = n+1;
                dir_matrix(:,n) = temp_dir(:,i);
            end
        end    
        
    case '10degOffset'
        for i = 1:Ndirs
            if (rem((temp_dir(1,i) + temp_dir(2,i)),10) == 0) & (rem(temp_dir(1,i),5)==0) & (rem(temp_dir(2,i),5)==0)
                n = n+1;
                dir_matrix(:,n) = temp_dir(:,i);
            end
        end
        
    case '10deg'
        for i = 1:Ndirs
            if (rem(temp_dir(1,i),10)==0) & (rem(temp_dir(2,i),10)==0)
                n = n+1;
                dir_matrix(:,n) = temp_dir(:,i);
            end
        end
        
    case '20degOffset'
        for i = 1:Ndirs
            if (rem((temp_dir(1,i) + temp_dir(2,i)-10),20) == 0) & (rem(temp_dir(1,i),10)==0) & (rem(temp_dir(2,i),10)==0)
                n = n+1;
                dir_matrix(:,n) = temp_dir(:,i);
            end
        end
        
    case '20deg'
        for i = 1:Ndirs
            if (rem(temp_dir(1,i),20)==0) & (rem(temp_dir(2,i)-10,20)==0)
                n = n+1;
                dir_matrix(:,n) = temp_dir(:,i);
            end
        end
        
    otherwise 
        disp ('ERROR: requested usupported spatial resolution')
end