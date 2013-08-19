function [filetype] = testHRTFfiletype(PATH, FN);

% function [filetype] = testHRTFfiletype(PATH, FN);
%
% function to test what type of HRTFfile this is
%   Filetype:
%       0) not yet tested
%       1) traditional binary
%       2) .mat (containing TF1 and TF2)
%       ...
%       999) bad file

dir = -1;

try
    eval(['load -mat ' PATH FN]);
    filetype = 2;
catch
    try
        eval(['dir = mtlrdir([PATH FN]);' ]);
        if dir ~=-1 
            filetype = 1;
        else
            disp([FN ' Unknown HRTF filetype or HRTF file does not exist'])
            filetype = 999;
        end
    catch
        disp([FN ' Unknown HRTF filetype or HRTF file does not exist'])
        filetype = 999;
    end
end