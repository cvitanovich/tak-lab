function [LT RT] = readHRTFdotMAT(fname,EL,AZ)
% [LT RT] = readHRTFdotMAT(fname,EL,AZ)
% grabs coefficients from an HRTF dot MAT file
load(fname)

% find index
[row col0]=find(dir(1,:)==EL);
[row col1]=find(dir(2,:)==AZ);
idx = intersect(col0,col1);

% get left and right coefficients
LT=TF1(idx,:);
RT=TF2(idx,:);

end

