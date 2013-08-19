function [LT,RT] = readHRTFdotMAT(hrtf_dir,fname,EL,AZ)
% [LT RT] = readHRTFdotMAT(fname,EL,AZ)
% grabs coefficients from an HRTF dot MAT file
load([hrtf_dir fname]);
% find index
EL_list=direc(1,:);
AZ_list=direc(2,:);
[row,col0]=find(EL_list==EL);
[row,col1]=find(AZ_list==AZ);
idx = intersect(col0,col1);

% get left and right coefficients
LT=TF1(idx,:);
RT=TF2(idx,:);

