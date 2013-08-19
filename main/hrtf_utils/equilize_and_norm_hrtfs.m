% Normalize and ABL Equalize an HRTF file
clear all
cd /Users/cvitanovich/Documents/MATLAB/tak-lab/HRTFs/Matlab_V6
uiopen('*.mat');

% ABL Equalize:
ABL_RMS=(sqrt(mean(TF1.^2,2))+sqrt(mean(TF2.^2,2)))/2;
mx=max(ABL_RMS);
for k=1:size(TF1,1)
    TF1(k,:)=(mx/ABL_RMS(k)).*TF1(k,:);
    TF2(k,:)=(mx/ABL_RMS(k)).*TF2(k,:);
end
% normalize so that output is always +/-1
absum_TF1=sum(abs(TF1),2);
absum_TF2=sum(abs(TF2),2);

mx2=max([absum_TF1; absum_TF2]);

TF1=(.9999*(1/mx2)).*TF1;
TF2=(.9999*(1/mx2)).*TF2;

HDR.comment3='ABL equalized and normalized for +/- 1 amplitudes.';
save('normed_hrtfs.mat','HDR','TF1','TF2','direc','-v6');
