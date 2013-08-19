function manual_fields(handle)
S=get(handle);
F=fieldnames(S);

for k=1:length(F)
    if(strcmp(S.(F{k}),'auto'))
        choices=set(handle,F{k});
        for j=1:length(choices)
            if(strcmp(choices{j},'manual'))
               set(handle,F{k},'manual');
            end
        end
    end
end