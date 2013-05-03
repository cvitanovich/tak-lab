function [trial_list] = randomized_trials(not_outlier_IDs,outlier_IDs,nOutliers,nTests,rep_limit)

outlier_list = [];
nOutlierTypes = length(outlier_IDs);
if nOutliers >0
    for i1=1:ceil(nOutliers/nOutlierTypes)
        outlier_list((i1-1)*nOutlierTypes+1:i1*nOutlierTypes)=outlier_IDs;
    end
    outlier_list = outlier_list(1:nOutliers);
end


not_outlier_list = [];
nNotOutlierTypes = length(not_outlier_IDs);
for i2=1:ceil((nTests-nOutliers)/nNotOutlierTypes)
    not_outlier_list((i2-1)*nNotOutlierTypes+1:i2*nNotOutlierTypes)=not_outlier_IDs;
end
not_outlier_list = not_outlier_list(1:(nTests-nOutliers));

total_list = [outlier_list not_outlier_list];

chk2=0;
cnt=0;
IDs = sort(unique(total_list));
while ~chk2
    chk2=1;
    trial_list = total_list(randperm(length(total_list)));
    for i2=outlier_IDs
        df_list = diff(find(trial_list==i2));
        if find(df_list < rep_limit)
            chk2=0;
        end
    end
    
    % check for repeated speakers
    for i2=IDs
        df_list = diff(find(trial_list==i2));
        if find(df_list == 1)
            chk2=0;
        end
    end
    cnt=cnt+1;
    if cnt>10000
        hWarn=warndlg('Cannot find a randomized trial list with these parameters!');
        uiwait(hWarn);
        trial_list=[];
        return;
    end
end