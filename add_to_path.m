% ADDS THIS CODE REPOSITORY TO THE MATLAB PATH
tmp=which('add_to_path');
tmp=tmp(1:end-13);
addpath(genpath(tmp));
% remove hidden git directories
p=path;
start=1;
flag_remove=0;
for j=1:length(p)
    if(j>=4)
        if(strcmp(p(j-3:j),'.git'))
            flag_remove=1;
        end
    end
    if(j>1 && strcmp(p(j-1),';'))
        start=j;
    end
    if(strcmp(p(j),';'))
        stop=j-1;
        if(flag_remove)
            rmpath(p(start:stop));
            flag_remove=0;
        end
    end
end