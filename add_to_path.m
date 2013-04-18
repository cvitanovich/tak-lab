% ADDS THIS CODE REPOSITORY TO THE MATLAB PATH
tmp=which('add_to_path');
tmp=tmp(1:end-13);
addpath(genpath(tmp));
% remove hidden git directory
rmpath(genpath([tmp '.git']));