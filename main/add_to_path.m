% ADDS THIS CODE REPOSITORY TO THE MATLAB PATH

current=pwd;
tmp=which('add_to_path');
tmp=tmp(1:end-13);
cd(tmp);
addpath(genpath(tmp));
cd(current);

% OLD don't use
if(0)
current=pwd;
tmp=which('add_to_path');
tmp=tmp(1:end-13);
cd(tmp);
addpath(genpath(tmp));
cd .git
rmpath(genpath(pwd));
cd(current);
end

% OLD APPROACH
if(0)
    % remove hidden git directories
    p=path;
    start=1;
    flag_remove=0;
    cmp_str=computer;
    if(strcmp(cmp_str(1:3),'MAC') || strcmp(cmp_str,'GLNXA64'))
        delimit=':';
    else
        delimit=';';
    end
end

if(0)
    for j=1:length(p)
        if(j>=4)
            if(strcmp(p(j-3:j),'.git'))
                flag_remove=1;
            end
        end
        if(j>1 && strcmp(p(j-1),delimit))
            start=j;
        end
        if(strcmp(p(j),delimit))
            stop=j-1;
            if(flag_remove)
                rmpath(p(start:stop));
                flag_remove=0;
            end
        end
    end
end